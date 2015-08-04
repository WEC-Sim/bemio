import subprocess
import os

class Version(object):

    def __init__(self):
        self._base = None
        self._full = None
        self.default = '1.1'

    # Get the code version from .git
    @ property
    def full(self, ):

        # try:
        #     ver = subprocess.Popen(["git", "describe","--tags", "--dirty", "--always"], stdout=subprocess.PIPE)
        #     self._full = ver.communicate()[0].rstrip()
        #
        # except:
        #     print "Unable to run Git on your system to determine the current bemio version"
        #     print 'Setting the version to the default: ' + self.default
        #     self._full = self.default

        self._full = self.default # hack
        return self._full

    @ property
    def base(self, ):

        self._base = self.full.split('-')[0]

        return self._base

def base():
    ver = Version()
    return ver.base

def full():
    ver = Version()
    return ver.full
