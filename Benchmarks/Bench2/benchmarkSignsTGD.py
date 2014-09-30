# -*- coding: iso-8859-15 -*-
"""
Benchmark Signs2 web application.


"""
import time
import unittest
from funkload.FunkLoadTestCase import FunkLoadTestCase
from webunit.utility import Upload


auto_refresh_string = 'This is an autorefreshing page'
MAX_running_time = 3600 * 1 


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
        self.server_url = 'http://signs2.iib.uam.es'
        ##self.server_url = self.conf_get('main', 'url')


    def breast(self):
        server_url = self.server_url

        self.get(server_url + "/",
            description="Get /")

        start_time = time.time()
        self.post(server_url + "/cgi-bin/signsR.cgi", params=[
            ['covariate', Upload("../DataSets/breast.covar.txt")],
            ['time', Upload("../DataSets/breast.surv.txt")],
            ['event', Upload("../DataSets/breast.event.txt")],
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


    def dlbcl(self):
        server_url = self.server_url

        self.get(server_url + "/",
            description="Get /")

        start_time = time.time()
        self.post(server_url + "/cgi-bin/signsR.cgi", params=[
            ['covariate', Upload("../DataSets/dlbcl.160.covar.txt")],
            ['time', Upload("../DataSets/dlbcl.160.surv.txt")],
            ['event', Upload("../DataSets/dlbcl.160.event.txt")],
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


    def breast_tgd(self):
        server_url = self.server_url

        self.get(server_url + "/",
            description="Get /")

        start_time = time.time()
        self.post(server_url + "/cgi-bin/signsR.cgi", params=[
            ['covariate', Upload("../DataSets/breast.covar.txt")],
            ['time', Upload("../DataSets/breast.surv.txt")],
            ['event', Upload("../DataSets/breast.event.txt")],
            ['methodSurv', 'TGD'],
            ['maxiter', '6000'],
            ['epi', '6e-06'],
            ['tau', '.90'],
            ['MinCor', '0.5'],
            ['organism', 'None'],
            ['idtype', 'None']],
            description="Post /cgi-bin/signsR.cgi")
        common_part_bench(self)

        end_time = time.time()
        duration = end_time - start_time
        print duration


    def dlbcl_tgd(self):
        server_url = self.server_url

        self.get(server_url + "/",
            description="Get /")

        start_time = time.time()
        self.post(server_url + "/cgi-bin/signsR.cgi", params=[
            ['covariate', Upload("../DataSets/dlbcl.160.covar.txt")],
            ['time', Upload("../DataSets/dlbcl.160.surv.txt")],
            ['event', Upload("../DataSets/dlbcl.160.event.txt")],
            ['methodSurv', 'TGD'],
            ['maxiter', '6000'],
            ['epi', '6e-06'],
            ['tau', '.90'],
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
