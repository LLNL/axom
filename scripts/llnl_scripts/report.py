#!/bin/sh
"exec" "python" "-u" "-B" "$0" "$@"

###############################################################################
# Copyright (c) 2017-2018, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-741217
#
# All rights reserved.
#
# This file is part of Axom.
#
# For details about use and distribution, please read axom/LICENSE.
###############################################################################

"""
 file: report.py

 description: 
  Emails an HTML report for jobs and test results

"""

import datetime
import operator
import os
import re
import sys

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET

from optparse import OptionParser
from smtplib import SMTP
from email.MIMEText import MIMEText
from dateutil.parser import parse
from llnl_lc_build_tools import *


__cellBlackOnGreen = "<td bgcolor=\"#00C000\"><font  color=\"#000000\">%s</font></td>"
__cellWhiteOnRed = "<td bgcolor=\"#E10000\"><font color=\"#FFFFFF\">%s</font></td>"
__contentHeader = "<html><head><meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\"></head><body>"
__contentFooter = "</body></html>"
__detailLinkIndex = 0


class JobInfo(object):
    def __init__(self):
        self.name = ""
        self.sys_type = ""
        self.archive_dir = ""
        self.success = False
        self.datetime = ""
        self.specInfos = {}


class SpecInfo(object):
    def __init__(self):
        self.name = ""
        self.success = False
        self.passed = []
        self.failed = []

def parse_args():
    "Parses args from command line"
    parser = OptionParser()
    # Location of source directory to build
    parser.add_option("-e", "--email",
                      type="string",
                      dest="email",
                      default="axom-dev@llnl.gov",
                      help="Email address to send report (Defaults to 'axom-dev@llnl.gov')")


    ###############
    # parse args
    ###############
    opts, extras = parser.parse_args()
    # we want a dict b/c the values could 
    # be passed without using optparse
    opts = vars(opts)
    return opts


def main():
    opts = parse_args()

    # Email info
    sender = "axom-dev@llnl.gov"
    receiver = opts["email"]
    emailServer = "nospam.llnl.gov"
    if on_rz():
        emailSubject = "[AXOM] RZ Status"
    else:
        emailSubject = "[AXOM] CZ Status"

    archive_dir = get_archive_base_dir()

    #Generate build email and send
    print "Reading archived job information..."
    basicJobInfos, specJobInfos, tplJobInfos = generateJobInfos(archive_dir)

    print "Sorting and filtering out old jobs..."
    basicJobInfos, specJobInfos, tplJobInfos = sortAndFilterJobInfos(basicJobInfos, specJobInfos, tplJobInfos)

    print "Generating email content..."
    emailContent = generateEmailContent(basicJobInfos, specJobInfos, tplJobInfos)
    
    print "Sending email..."
    return sendEmail(emailContent, emailSubject, sender, receiver, emailServer)


def getAllDirectoryNames(path):
    return [n for n in os.listdir(path) if os.path.isdir(pjoin(path, n))]


def determineSuccessState(path):
    success = False
    if os.path.exists(pjoin(path, "success.json")):
        success = True
    elif os.path.exists(pjoin(path, "failed.json")):
        success = False
    else:
        success = None
    return success


def generateJobInfos(archive_dir):
    basicJobInfos = {} # key's are job names, value is a list of jobs
    specJobInfos = {} # key's are sys_type, value is a list of jobs
    tplJobInfos = {} # key's are sys_type, value is a list of jobs
    sys_types = getAllDirectoryNames(archive_dir)

    # Directory structure = <archive base>/<sys_type>/<job name>/<datetime>/<spec>/[Logs]
    for sys_type in sys_types:
        sys_type_dir = pjoin(archive_dir, sys_type)
        jobNames = getAllDirectoryNames(sys_type_dir)
        for jobName in jobNames:
            jobName_dir = pjoin(sys_type_dir, jobName)
            datetimes = getAllDirectoryNames(jobName_dir)

            for datetime in datetimes:
                currJobInfo = JobInfo()
                currJobInfo.name = jobName
                currJobInfo.sys_type = sys_type
                currJobInfo.datetime = datetime
                datetime_dir = pjoin(jobName_dir, datetime)
                currJobInfo.archive_dir = datetime_dir

                isTPLJob = False
                specs = getAllDirectoryNames(datetime_dir)
                currJobInfo.success = determineSuccessState(datetime_dir)

                # Basic jobs should contain no folders
                if len(specs) == 0:
                    # Used for general purpose jobs
                    jobType = "basic"
                else:
                    # TPL jobs will have spec dirs and an info.json at the root
                    if os.path.exists(pjoin(datetime_dir, "info.json")):
                        jobType = "tpl"
                    else:
                        jobType = "spec"

                    # spec or tpl job
                    for spec in specs:
                        currSpecInfo = SpecInfo()
                        currSpecInfo.name = spec
                        spec_dir = pjoin(datetime_dir, spec)
                        currSpecInfo.success = determineSuccessState(spec_dir)

                        populateTests(currSpecInfo, spec_dir)
                        if currJobInfo.specInfos.has_key(currSpecInfo.name):
                            print "Warning: duplicate spec ({0}) found in job: {1}".format(spec, currJobInfo.name)
                        currJobInfo.specInfos[spec] = currSpecInfo

                        currJobInfo.success = True
                        for specName in currJobInfo.specInfos.keys():
                            # if any tests have failed then the whole job fails
                            currSpecInfo = currJobInfo.specInfos[specName]
                            if len(currSpecInfo.failed) > 0:
                                currJobInfo.success = False
                                break

                if jobType == "basic":
                    if not basicJobInfos.has_key(jobName):
                        basicJobInfos[jobName] = []
                    basicJobInfos[jobName].append(currJobInfo)
                elif jobType == "spec":
                    if not specJobInfos.has_key(currJobInfo.sys_type):
                        specJobInfos[currJobInfo.sys_type] = []
                    specJobInfos[currJobInfo.sys_type].append(currJobInfo)
                elif jobType == "tpl":
                    if not tplJobInfos.has_key(currJobInfo.sys_type):
                        tplJobInfos[currJobInfo.sys_type] = []
                    tplJobInfos[currJobInfo.sys_type].append(currJobInfo)
                else:
                    print "Error: Unsupported job type: {0}".format(jobType)

    return basicJobInfos, specJobInfos, tplJobInfos


def populateTests(specInfo, path):
    test_xml_path = pjoin(path, 'Test.xml')
    if not os.path.exists(test_xml_path):
        return

    #xml_root = ET.parse(test_xml_path).getroot()
    tree = ET.parse(test_xml_path)

    for test_elem in tree.iter(tag="Test"):
        # Skip the list of test names @ Site -> Testing -> TestList -> Test
        if not test_elem.attrib.has_key("Status"):
            continue

        name = ""
        for child_elem in test_elem:
            if child_elem.tag == "Name":
                name = child_elem.text
                break
        if name == "":
            print "Error: {0}: Unknown to find test name".format(test_xml_path)

        if test_elem.attrib["Status"] == "passed":
            specInfo.passed.append(name)
        elif test_elem.attrib["Status"] == "failed":
            specInfo.failed.append(name)
        else:
            print "Error: {0}: Unknown test status ({1})".format(test_xml_path, test_elem.attrib["Status"])


def sortAndFilterJobInfos(basicJobInfos, specJobInfos, tplJobInfos):
    # Sort all jobs with the same name with newest job first based on date time
    # and cut off at the number we want to report
    for jobName in basicJobInfos.keys():
        # we just want the newest one
        basicJobInfos[jobName].sort(key=operator.attrgetter("datetime"))
        basicJobInfos[jobName].reverse()
        del basicJobInfos[jobName][1:]

    for sys_type in specJobInfos.keys():
        # we just want the newest 5
        jobInfoList = specJobInfos[sys_type] # This gets all jobs with the specific sys_type
        jobInfoList.sort(key=operator.attrgetter("datetime"))
        jobInfoList.reverse()
        del jobInfoList[5:]

    for sys_type in tplJobInfos.keys():
        # we just want the newest 5
        jobInfoList = specJobInfos[sys_type] # This gets all jobs with the specific sys_type
        jobInfoList.sort(key=operator.attrgetter("datetime"))
        jobInfoList.reverse()
        del jobInfoList[5:]

    return basicJobInfos, specJobInfos, tplJobInfos


def generateEmailContent(basicJobInfos, specJobInfos, tplJobInfos):
    html = __contentHeader

    # Add Basic jobs to the email
    if len(basicJobInfos.keys()) > 0:
        # Header for all basic jobs
        html += "<center><font size=\"5\"> </font></center><br>"
        html += "<center><b><font size=\"5\">Basic Jobs</font></b></center>"
        html += "<table border=\"1\">"

        # Start table for basic jobs
        html += "<tr>"
        html += "<th bgcolor=\"#C0C0C0\">Name</th>"
        html += "<th bgcolor=\"#C0C0C0\">Date</th>"
        html += "<th bgcolor=\"#C0C0C0\">Success</th>"
        html += "</tr>\n"

        names = basicJobInfos.keys()
        names.sort()
        for name in names:
            currJobInfo = basicJobInfos[name][0]
            # Row for each job
            html += "<tr>"
            html += "<td>" + currJobInfo.name + "</td>"
            html += "<td>" + convertDatetimeToReadable(currJobInfo.datetime) + "</td>"
            html += getSuccessFormat(currJobInfo.success, __cellBlackOnGreen, __cellWhiteOnRed) % str(currJobInfo.success)
            html += "</tr>"

        # Close table for basic jobs
        html += "</table>\n"

    html = getHTMLforJobInfos(specJobInfos, html, False)
    html = getHTMLforJobInfos(tplJobInfos, html, True)

    html += __contentFooter
    return html


def getHTMLforJobInfos(jobInfosDict, html, isTPLJob):
    if len(jobInfosDict.keys()) > 0:
        # Header for all spec jobs
        html += "<center><font size=\"5\"> </font></center><br>"
        if (isTPLJob):
            html += "<center><b><font size=\"5\">TPL Jobs</font></b></center>"
        else:    
            html += "<center><b><font size=\"5\">Spec Jobs</font></b></center>"

        sys_types = jobInfosDict.keys()
        sys_types.sort()
        for sys_type in sys_types:
            jobInfos = jobInfosDict[sys_type]
            specNames = getSpecNames(jobInfos)

            # Header for this sys_type
            html += "<center><font size=\"5\"> </font></center><br>"
            html += "<font size=\"3\">" + sys_type + "</font>"

            # Start table for this spec jobs
            html += "<table border=\"1\">"
            html += "<tr>"
            html += "<th bgcolor=\"#C0C0C0\">Date</th>"
            html += "<th bgcolor=\"#C0C0C0\">Success</th>"
            for specName in specNames:
                html += "<th bgcolor=\"#C0C0C0\">" + specName + "</th>"
            html += "</tr>\n"

            for currJobInfo in jobInfos:
                # Row for each job
                failedTestsString = ""
                html += "<tr>"
                html += "<td>" + convertDatetimeToReadable(currJobInfo.datetime) + "</td>"
                html += getSuccessFormat(currJobInfo.success, __cellBlackOnGreen, __cellWhiteOnRed) % str(currJobInfo.success)
                for specName in specNames:
                    if currJobInfo.specInfos.has_key(specName):
                        currSpecInfo = currJobInfo.specInfos[specName]
                        totalNumberOfTests = len(currSpecInfo.passed) + len(currSpecInfo.failed)
                        cellString = "{0}/{1}".format(str(len(currSpecInfo.passed)), str(totalNumberOfTests))
                        html += getSuccessFormat(currSpecInfo.success, __cellBlackOnGreen, __cellWhiteOnRed) % cellString

                        # Build up a list of failed tests to be added to the html table at the end
                        cnt = 0
                        if len(currSpecInfo.failed) > 0:
                            if failedTestsString != "":
                                failedTestsString += "<br />"
                            failedTestsString += "[{0}]".format(currSpecInfo.name)
                            for name in currSpecInfo.failed:
                                if failedTestsString[-1] != "]":
                                    failedTestsString += ",\n" if cnt % 20 == 0 else ","
                                failedTestsString += " " + name
                                cnt += 1
                    else:
                        # Have a blank cell if spec wasn't here for this job (it was added or removed)
                        html += "<td></td>"
                html += "</tr>"
                tableWidth = 2 + len(specNames)
                if failedTestsString != "":
                    html += "<tr><td colspan=\"{0}\">{1}</td></tr>".format(str(tableWidth), failedTestsString)
                html += "<tr><td colspan=\"{0}\">{1}</td></tr>".format(str(tableWidth),currJobInfo.archive_dir)

            # Close table for this sys_type
            html += "</table>\n"

    return html


def getSpecNames(jobInfos):
    specNames = []

    for currJobInfo in jobInfos:
        for specName in currJobInfo.specInfos.keys():
            specNames.append(specName)
    specNames = list(set(specNames)) # keep only unique names
    specNames.sort()

    return specNames


def convertDatetimeToReadable(datetimestring):
    import datetime
    dt = datetime.datetime.strptime(datetimestring, "%Y_%m_%d_%H_%M_%S")
    dateString = dt.strftime("%a, %d %b %Y %H:%M:%S")
    return dateString


def getSuccessFormat(success, successFormat, failureFormat):
    if success:
        return successFormat
    return failureFormat


def sendEmail(content, subject, sender, reciever, emailServer):
    try:
        msg = MIMEText(content, "html")
        msg["Subject"] = subject
        msg["From"] = sender
        msg["To"] = reciever
        msg["Content-Type"] = "text/html"

        conn = SMTP(emailServer)
        conn.set_debuglevel(False)
        try:
            conn.sendmail(sender, reciever, msg.as_string())
        finally:
            conn.close()
    except Exception as e:
        print "Failed to send email:\n {0}".format(str(e))
        return False
    print "Sent."

    return True


if __name__ == "__main__":
    if main():
        sys.exit(0)
    sys.exit(1)

