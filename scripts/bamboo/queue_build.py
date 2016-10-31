#! /usr/bin/env python

#  Author:  Ian Lee, 28 June 2016
#  https://mystash.llnl.gov/projects/ATLASSIAN/repos/scripts/browse/queue_build.py

import getpass
import logging
import requests
import getopt, sys

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def main(argv):
   global testPlan
   try:
      opts, args = getopt.getopt(argv,"h:p",["plan="])
   except getopt.GetoptError:
      print 'queue_build.py  -p <plan> '
      sys.exit(2)
#   print "in main"
   for opt, arg in opts:
      if opt == '-h':
         print 'queue_build.py  -p <plan> '
         sys.exit()
      elif opt in ("-p", "--plan"):
         testPlan = arg
         print testPlan
      else:
         print 'queue_build.py  -p <plan> '
         sys.exit()


__url_cache__ = {}

class CZBamboo(requests.Session):

    def __init__(self):

        super(CZBamboo, self).__init__()
        self.headers.update({
            # Only accept JSON responses
            'Accept': 'application/json',
            # Only accept UTF-8 encoded data
            'Accept-Charset': 'utf-8',
            # Always sending JSON
            'Content-Type': "application/json",
        })

        self.login_url = 'https://lc.llnl.gov/dologin.cgi'
        self.base_url = 'https://lc.llnl.gov/bamboo/rest/api/latest'

    def build_url(self, *args, **kwargs):
        """
        Builds a new API url from scratch.
        Adapted from: https://github.com/sigmavirus24/github3.py/blob/develop/github3/session.py
        """

        parts = [kwargs.get('base_url') or self.base_url]
        parts.extend(args)
        parts = [str(p) for p in parts]
        key = tuple(parts)
        logger.info('Building a url from %s', key)

        if key not in __url_cache__:
            logger.info('Missed the cache building the url')
            __url_cache__[key] = '/'.join(parts)
        return __url_cache__[key]

    def login(self):
        """
        Login to Lorenz with credentials gained via terminal prompt
        """

        username = getpass.getuser()
        password = getpass.getpass('Pin & Token: ')
        response = self.post(self.login_url, auth=(username, password))
        logger.debug('Server response: %s', response.__dict__)

    def queue(self, plan_key=None):

        if plan_key is not None:
            self.post(self.build_url('queue', plan_key))
        response = self.get(self.build_url('queue'))
        return response.json()


if __name__ == '__main__':

    main(sys.argv[1:])
    msg = "Test plan = " + testPlan
#    print msg
    bamboo = CZBamboo()
    bamboo.login()
#    queue = bamboo.queue('ASC-NIG')
#    queue = bamboo.queue('ASC-NIG71')
    queue = bamboo.queue(testPlan)
    print(queue)


