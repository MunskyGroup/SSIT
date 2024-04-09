#!/bin/bash
# Author: Brian Munsky
# Date: 09/24/16
# Purpose: Batch Matlab Job Submission

import subprocess
import time
for i in range(1,5,1): 
	for j in range(1,3,1):  
		cmd = ' '.join( ['qsub','-q musky.q@node* -cwd -o /dev/null -e Errs/','-v','','i=%d' % (i),'-v','','j=%d' % (j),'Python_Wrapper.sh'] )
		print cmd
		subprocess.call( cmd, shell=True )
		time.sleep(0.1) # delay for 0.1 seconds.	
