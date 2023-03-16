# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:07:09 2023

@author: gegese
"""

import requests


def brute_downloader(basepath):
    
    all_ids = []
    
    printcounter = 0
    
    for i in range(200000,400000):
        
        if (printcounter == 500):
            print('At iteration %s' % i)
            printcounter = 0
        
        r = requests.get(r'https://ehyd.gv.at/eHYD/MessstellenExtraData/owf?id=%s&file=1' % i, allow_redirects=True)
        r2 = requests.get(r'https://ehyd.gv.at/eHYD/MessstellenExtraData/owf?id=%s&file=4' % i, allow_redirects=True)
        
        if len(r.content) > 0.:
            open(basepath + '/%s_meta.csv' % i, 'wb').write(r.content)
            open(basepath + '/%s_runoff.csv' % i, 'wb').write(r2.content)
    
    return all_ids
    
basepath = r'C:\Projects\DataMining\ZAMG_API\Gauges'

brute_downloader(basepath)
