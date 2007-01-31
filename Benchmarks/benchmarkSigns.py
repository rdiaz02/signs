# -*- coding: iso-8859-15 -*-
"""
Benchmark SignS web application.


"""
import time
import unittest
from funkload.FunkLoadTestCase import FunkLoadTestCase
from webunit.utility import Upload


auto_refresh_string = 'This is an autorefreshing page'
MAX_running_time = 3600 * 1 


def common_part(self, final_output,                
                MAX_running_time = 3600,
                auto_refresh_string = auto_refresh_string):    
    server_url = self.server_url
    start_run = time.time()
    refresh_num = 0
    
    while True:
        final_body = self.getBody()
        if final_body.find(auto_refresh_string) < 0:
            break
        time.sleep(5)
        refresh_num += 1
        run_time = time.time() - start_run
        print '\n Refreshed ' + str(refresh_num) + ' times. Been running for ' + str(round(run_time/60.0, 2)) + ' minutes.\n'
        if run_time > MAX_running_time :
            self.fail('Run longer than MAX_running_time')
        self.get(server_url + self.getLastUrl(),
                 description="Get /cgi-bin/checkdone.cgi")
    expected = final_body.find(final_output) >= 0
    if not expected:
        self.fail('\n ***** (begin of) Unexpected final result!!!! *****\n' + \
                 str(final_body) + \
                 '\n ***** (end of) Unexpected final result!!!! *****\n')
    else:
        print 'OK'


def common_part_bench(self,                 
                      MAX_running_time = 3600,
                      auto_refresh_string = auto_refresh_string):    
    """ like above, but does not check anything. simply benchmarking"""
    server_url = self.server_url
    
    while True:
        final_body = self.getBody()
        if final_body.find(auto_refresh_string) < 0:
            break
        time.sleep(5)
        self.get(server_url + self.getLastUrl(),
                 description="Get /cgi-bin/checkdone.cgi")
    print 'OK'

   
class Signs(FunkLoadTestCase):
    """XXX

    This test use a configuration file Signs.conf.
    """

    def setUp(self):
        """Setting up test."""
        self.logd("setUp")
        self.server_url = 'http://signs.bioinfo.cnio.es'
        ##self.server_url = self.conf_get('main', 'url')

    def test1(self):
        server_url = self.server_url

        self.get(server_url + "/",
            description="Get /")

        self.post(server_url + "/cgi-bin/signsR.cgi", params=[
            ['covariate', Upload("./mini.lung.covar.txt")],
            ['time', Upload("./mini.lung.surv.txt")],
            ['event', Upload("./mini.lung.event.txt")],
            ['methodSurv', 'FCMS'],
            ['Minp', '0.1'],
            ['MaxSize', '100'],
            ['MinSize', '10'],
            ['MinCor', '0.5'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="Post /cgi-bin/signsR.cgi")
        final_output = 'No groups that meet the p, minimum correlation and size restrictions.'
        common_part(self, final_output)


    def test1b(self):
        server_url = self.server_url

        self.get(server_url + "/",
            description="Get /")

        start_time = time.time()
        self.post(server_url + "/cgi-bin/signsR.cgi", params=[
            ['covariate', Upload("./mini.lung.covar.txt")],
            ['time', Upload("./mini.lung.surv.txt")],
            ['event', Upload("./mini.lung.event.txt")],
            ['methodSurv', 'FCMS'],
            ['Minp', '0.1'],
            ['MaxSize', '100'],
            ['MinSize', '10'],
            ['MinCor', '0.5'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="Post /cgi-bin/signsR.cgi")
        common_part_bench(self)

        end_time = time.time()
        duration = end_time - start_time
        print duration

    def tearDown(self):
        """Setting up test."""
        self.logd("tearDown.\n")

if __name__ in ('main', '__main__'):
     unittest.main()
