import os, sys
from commonOps import *

# Class for running jobs on the cluster
class ClusterJobs(object):
    def __init__(self, config, runDir, templateStr, list1, list2):
        # Name of cluster to use, config is the parsed Config file 
        self.cluster = config.getGenVar('cluster')
        # Template file for creating list of jobs
        self.template = self.__createTemplate(runDir, templateStr)
        # File with list of jobs to run on the cluster
        self.jobList = createPath(runDir, "jobList")
        # List of items to be substituted into template, full path
        self.list1 = list1
        # Second list of items to be substituted into template, full path.
        # This is "single" if there is only one list.
        self.list2 = list2
        # Program and list of arguments to run jobs
        batchOption = "-batch=" + runDir
        self.makeJobs = ["para", "make", self.jobList, batchOption, "-ram=4g"] 

    def __createTemplate(self, runDir, templateStr):
        """Creates a file with the template for substitution of file lists
           to create a list of jobs for the cluster run"""
        templateFile = createPath(runDir, "template")
        with open (templateFile, "w") as tFh:
            tFh.write("#LOOP\n")
            tFh.write(templateStr + "\n")
            tFh.write("#ENDLOOP\n")
        # Return the template file name
        return templateFile

    def createJobList(self):
        """Substitutes items from list(s) into template to create jobList"""
        subprocess.check_call(["gensub2", self.list1, self.list2, \
            self.template, self.jobList])

    def runJobs(self):
        """Creates list of jobs and runs the jobs in the jobList"""
        self.createJobList()
        # Run jobs on cluster
        subprocess.check_call(["ssh", "-T", self.cluster] + self.makeJobs) 
        # Creates a job database and manages running jobs on cluster
        subprocess.check_call(self.makeJobs)

