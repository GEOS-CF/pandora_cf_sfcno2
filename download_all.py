#!/usr/local/other/python/GEOSpyD/2019.03_py3.7/2019-04-22/bin/python
'''
Helper script to loop through all entries of the provided json file
and attempt to download the file provided from key 'pandora_url'.
'''
import os
import argparse
import json
import urllib.request

def main(args):
    urls = json.load(open(args.urllist))
    for l in urls:
        url = l.get('pandora_url',None)
        if url is None:
            continue
        ofile = url.split('/')[-1]
        if not os.path.isfile(ofile):
            print('attempting to download {}'.format(url))
            try:
                urllib.request.urlretrieve(url,ofile)
            except:
                print('could not download {}'.format(ofile))
        else:
            print('file exists - skip: {}'.format(ofile))
    return


def parse_args():
    p = argparse.ArgumentParser(description='Undef certain variables')
    p.add_argument('-l', '--urllist',type=str,help='url list',default="PANDORA_Locations.json")
    return p.parse_args()

 
if __name__ == '__main__':
    main(parse_args())
#
