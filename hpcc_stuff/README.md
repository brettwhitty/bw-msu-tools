## HPCC scripts

### exonerate\_submitter.pl

This is a crazy workaround for the fact the gateway host for the MSU campus HPCC was not a grid submit host.

Consequence of this was that users were required to SSH to gateway, then SSH to the submit host to submit jobs.

I wrote this script to be able to do command-line job submissions from the buell-lab server by having this script copy itself across hosts before doing a job submission.

### pbs\_tools/

Experiment building OO tool modules on top of Perl 'PBS::Client' module. See README.md in directory.
