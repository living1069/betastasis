#!/tools/bin/python

# NIH TCGA downloader

import httplib2
import os
import sys
import time
import errno
import re
import ConfigParser
import threading
import Queue

queue = Queue.Queue()

dir_skip_patterns = [
        #'^minbio',
        #'^bio$',
        #'^tracerel$',
        'tissue_images',
        #'.*Level_2.*',
        #'.*Level_3.*',
        #'.*\.aux\..*',
        #'archive',
        #'hg-cgh-415k_g4124a',
        #'mskcc.org',
        #'^gsc$',
        #'humanhap',

        # These are redundant old-style archives whose contents are also available
        # in unified Level_1 archives.
        'unc.edu_GBM.AgilentG4502A_07_2.\d',
        'unc.edu_OV.AgilentG4502A_07_3.\d',

        # These rules were used to limit downloads to only Agilent gene expression
        # data.
        #'hms.harvard.edu',
        #'broad.mit.edu',
        #'jhu-usc.edu',
        #'mskcc.org',
        #'hudsonalpha.org',
]

file_skip_patterns = [
        #'.*level2.*',
        #'.*level3.*',
        '.*\.zip$',
        '.*\.tar\.gz$',
        '.*\.jpg$',
        '.*\.JPG$',
        '.*\.png$',
        '.*\.pdf$',
        '.*\.md5$',
        '.*\.tif$',
        '.*\.TIF$',
        #'.*\.tsv$',
        '.*\.gif$',
        #'.*_genes\.txt$',
        #'.*_tags\.txt$',
        #'.*_frequency\.txt$',
        'nohup.out',
]

class DownloadInfo(object):
    def __init__(self, url, path):
        self.url = url
        self.path = path

class DownloadThread(threading.Thread):
    def __init__(self, **kwargs):
        threading.Thread.__init__(self)

        self.queue = queue
        self.conn = httplib2.Http()

        if 'username' in kwargs and 'password' in kwargs:
            self.conn.add_credentials(kwargs["username"], kwargs["password"])

    def run(self):
        while True:
            dwninfo = self.queue.get()

            resp, content = self.conn.request(dwninfo.url, 'GET')
            if resp['status'] != '200':
                print resp
                print 'HTTP request failed for ' + dwninfo.url
                return

            listing = []
            for line in content.splitlines():
                m = re.match('^\s*<a href=.+?>(.+)</a>.*', line)
                if not m: continue
                listing.append(m.group(1))

            subdirs = []
            files = []

            for filename in listing:
                if filename[-1] == '/':
                    subdirs.append(filename[0:-1])
                else:
                    files.append(filename)

            for ffile in filter_files(files):
                ftp_size = file_size(self.conn, dwninfo.url + ffile)
                ffile_path = dwninfo.path + '/' + ffile

                local_size = 0
                if not os.path.exists(ffile_path):
                    if os.path.exists(ffile_path + '.bz2'):
                        local_size = ftp_size
                    elif os.path.exists(ffile_path + '.gz'):
                        local_size = ftp_size
                else:
                    local_size = os.path.getsize(ffile_path)

                if local_size != 0 and local_size == ftp_size:
                    print 'file exists %s' % (dwninfo.path + '/' + ffile)
                    sys.stdout.flush()
		    time.sleep(1)
                    continue

                print '+ url %s time[%s]' % (dwninfo.url + ffile, time.strftime("%c"))
                sys.stdout.flush()
		time.sleep(1)
		resp, content = self.conn.request(dwninfo.url + '/' +  ffile, 'GET')
            	if resp['status'] != '200':
                	print resp
                	print 'HTTP GET failed for ' + dwninfo.url + ffile
                	return

                localFile = open('%s/%s' % (dwninfo.path, ffile), 'w')
                localFile.write(content)
                localFile.close()
		print 'wrote %s at time[%s]' % (dwninfo.path + ffile, time.strftime("%c"))
            for fsub in filter_subdirs(subdirs):
                handle_subdir(dwninfo.url + '/' + fsub, dwninfo.path + '/' + fsub)

            self.queue.task_done()
	    #return

def filter_files(files):
    filtered = []
    for file in files:
        skip = False
        for pat in file_skip_patterns:
            if re.match(pat, file):
                skip = True
                break
        if not skip: filtered.append(file)
    return filtered

def filter_subdirs(dirs):
    pre_filtered = []
    filtered = []

    for dir in dirs:
        skip = False
        for pat in dir_skip_patterns:
            if re.match(pat, dir):
                skip = True
                break
        if not skip: pre_filtered.append(dir)

    versioned = []
    prefix = []
    minor = []

    for dir in pre_filtered:
        m = re.match('(.+)\.Level_(1|2|3)\.(\d+)\.(\d+)\.\d+$', dir)
        if m:
            versioned.append(dir)
            prefix.append(m.group(1) + '.' + m.group(2))
            minor.append(int(m.group(3)))
            continue

        m = re.match('(.+\.\d+)\.(\d+)\.\d+$', dir)
        if m:
            versioned.append(dir)
            prefix.append(m.group(1))
            minor.append(int(m.group(2)))
            continue

        filtered.append(dir)

    for m in set(prefix):
        highest_minor = -1
        highest_minor_idx = -1
        for n in range(len(prefix)):
            if prefix[n] != m: continue
            if minor[n] > highest_minor:
                highest_minor = minor[n]
                highest_minor_idx = n

        filtered.append(versioned[highest_minor_idx])

    return filtered

def file_size(conn, url):
    try:
    	resp, content = conn.request(url, 'HEAD', headers={'Accept-Encoding': 'plain'})
    	return int(resp['content-length'])
    except KeyError:
	#return some random number if content length not found
    	return 9999

def handle_subdir(targethost, localpath):
    targethost.rstrip('/')
    localpath.rstrip('/')

    time.sleep(1)
    print 'Processing dir: %s' % (localpath)

    try:
        os.makedirs(localpath)
    except OSError, exc:
        if exc.errno == errno.EEXIST: pass
        else: raise

    queue.put(DownloadInfo(targethost, localpath))

    return

def main(argv):
    if len(sys.argv) != 3:
        print 'Illegal arguments. Example of proper usage: tcga_download.py <config file> <tumor type>'
        sys.exit(-1)

    config = ConfigParser.RawConfigParser()
    config.read(sys.argv[1])

    subFilePath = sys.argv[2]

    host = config.get("Connection", "host")
    host.rstrip('/')

    username = config.get("Connection", "username")
    password = config.get("Connection", "password")

    localPath = config.get("Local", "localPath")
    localPath.rstrip('/')

    numberOfThreads = config.getint("Local", "numberOfThreads")

    print 'Mirroring data for tumor type %s start [%s]' % (subFilePath, time.strftime("%c"))
    print ('  from %s' % host)
    print ('    to %s' % localPath)

    print 'Starting Thread Pool with %s threads' % numberOfThreads
    for i in range(numberOfThreads):
        ts = DownloadThread(username=username, password=password)
        ts.setDaemon(True)
        ts.start()

    handle_subdir(host + '/' + subFilePath, localPath + '/' + subFilePath)
    print 'Mirroring data for tumor type %s completed [%s]' % (subFilePath, time.strftime("%c"))

    queue.join()

if __name__ == "__main__":
    main(sys.argv[1:])

