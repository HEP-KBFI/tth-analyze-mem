#!/usr/bin/env python
import subprocess, sys, getpass, time

if __name__ == '__main__':
  commands, sbatchIDs = [], []
  {% for bashScript, logFile in zippedScriptLog %}
  commands.append('sbatch --output={{ logFile }} {{ bashScript }}') {% endfor %}
  for command in commands:
    submitJobProcess = subprocess.Popen(command,
                                        stdout = subprocess.PIPE,
                                        stderr = subprocess.PIPE,
                                        shell = True)
    submitStdout, submitStderr = submitJobProcess.communicate()
    sbatchIDs.append(submitStdout.rstrip('\n').split()[-1])

  squeueCommand = "squeue -u %s | tail -n+2 | awk '{print $1}'" % getpass.getuser()
  while True:
    squeueProcess = subprocess.Popen(squeueCommand,
                                     stdout = subprocess.PIPE,
                                     stderr = subprocess.PIPE,
                                     shell = True)
    squeueStdout, squeueStderr = squeueProcess.communicate()
    nofJobsStillRunning = 0
    for line in squeueStdout.rstrip('\n').split('\n'):
      if line in sbatchIDs: nofJobsStillRunning += 1
    if nofJobsStillRunning == 0: break
    else:                        time.sleep(30)

