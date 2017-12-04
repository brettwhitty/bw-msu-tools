## pbs\_tools

This code was written as a test of using the Perl PBS::Client module to develop my own PBS::Tool::* modules for running multi-step bioinformatics tools on the MSU HPCC (via PBS Torque) in a simple way.

I created 'PBS::Tool::AAT' as a test case, but was dissatified with the real-world performance and did not end up using in production.

See:
	
	./bin/pbs_tool_test.pl

for an example use case of the AAT tool module.