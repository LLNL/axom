#
# Query Python executable
# Return semicolon delimited list of values.
#   version ; prefix
#
import sys

output = []
output.append(str(sys.version_info.major))
output.append(str(sys.version_info.minor))
output.append(str(sys.version_info.micro))
output.append(sys.prefix)

#sys.stdout.write('"' + ';'.join(output) + '"')
sys.stdout.write(';'.join(output) + "\n")
