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

queue = Queue.Queue(0)
fileQueue = Queue.Queue(1)
errorLogFile = None
numberOfDirThreads = 0
numberOfFileThreads = 0
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

dir_include_patterns = ['.*Level_*', '.*aux*', '.*mage-tab*']

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


class DownloadFile(object):
    def __init__(self, url, file):
       self.url = url
       self.file = file

class DownloadFileThread(threading.Thread):

	def __init__(self, **kwargs):
		threading.Thread.__init__(self)
		self.conn = httplib2.Http()

		if 'username' in kwargs and 'password' in kwargs:
			self.conn.add_credentials(kwargs["username"], kwargs["password"])

	def run(self):
		#check is isEmpty
		#while True and (not fileQueue.empty()):
		while True:    
			fileinfo = fileQueue.get()
			if fileinfo != None: 
				time.sleep(1)
				resp, content = self.conn.request(fileinfo.url)
				if resp['status'] == '200':
					localFile = open('%s' % (fileinfo.file), 'w')
					localFile.write(content)
					localFile.close()
					print 'Downloaded %s: \n%s' % (time.strftime("%c"), fileinfo.file)
				else:
					print 'DownloadFileThread: Error Getting %s  \n Resp: %s' % (fileinfo.url, resp)
					errorLogFile.write("Error Downloading url " + fileinfo.url + '\n')
			#for i in range(2):
			#    threadName = 'Dir' + str(i)
			#    if thread.Thread(threadName).isAlive():	       	         
			fileQueue.task_done()
			#if fileQueue.empty() and dirQueue.empty():
			#	print 'Queues are empty, downloading job completed %s' % (time.strftime("%c")) 
     			
  

class DownloadDir(object):
    def __init__(self, url, path):
        self.url = url
        self.path = path


class DownloadDirThread(threading.Thread):
    def __init__(self, **kwargs):
        threading.Thread.__init__(self)

        self.queue = queue
        self.conn = httplib2.Http()

        if 'username' in kwargs and 'password' in kwargs:
	    #print 'username pw %s %s' % (kwargs["username"], kwargs["password"])		
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
                ffile_path = dwninfo.path + '/' + ffile
                local_size = 0
                #if file is not in local path, add to queue
		if not os.path.exists(ffile_path):
		    fileQueue.join()	
		    fileQueue.put(DownloadFile(dwninfo.url + '/' + ffile, ffile_path))
		    print '++ new resource to fileQueue %s' % (dwninfo.url + '/' + ffile)
		    continue		
                    #if os.path.exists(ffile_path + '.bz2'):
                    #    local_size = ftp_size
                    #elif os.path.exists(ffile_path + '.gz'):
                    #    local_size = ftp_size
                else:
                    local_size = os.path.getsize(ffile_path)

		ftp_size = file_size(self.conn, dwninfo.url + '/' + ffile)
                #ffile_path = dwninfo.path + '/' + ffile
                #local_size = 0
                if local_size != 0 and local_size == ftp_size:
                	print 'Skip file %s local version has same file size %s' % (ffile_path, str(local_size))
                	sys.stdout.flush()
			#time.sleep(1)			
                	#continue
		else:
			fileQueue.join()
			fileQueue.put(DownloadFile(dwninfo.url + '/' + ffile, ffile_path))
                        print '++ update to fileQueue %s' % (dwninfo.url + '/' + ffile)

            for fsub in filter_subdirs(subdirs):
                handle_subdir(dwninfo.url + '/' + fsub, dwninfo.path + '/' + fsub)

            self.queue.task_done()
	    #return
	    #time.sleep(1)	
	    #if fileQueue.empty() and self.queue.empty():
            #	print 'Queues are empty, downloading job completed %s' % (time.strftime("%c"))
	
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

    for dir in dirs:
        skip = False
        for pat in dir_skip_patterns:
            if re.match(pat, dir):
                skip = True
                break
        if not skip: pre_filtered.append(dir)

    return pre_filtered


def file_size(conn, url):
	try:	
		resp, content = conn.request(url, 'HEAD', headers={'Accept-Encoding': 'plain'})
		return int(resp['content-length'])
	except KeyError:
		print "Key error on checking remote file content-length " + resp
		errorLogFile.write("Error getting file size for " + url + '\n')
		return -1

def make_path(localpath):
	try:
		os.makedirs(localpath)
	except OSError, exc:
		if exc.errno == errno.EEXIST: pass
		else: raise


def handle_subdir(targethost, localpath):
	targetUrlLevel = targethost.split('/')
	level = len(targetUrlLevel)
	targethost.rstrip('/')
	localpath.rstrip('/')

	#time.sleep(1)
	print 'Process dir: %s level %s' % (targethost, str(level))
	fileQueue.join()
	#download files from only directories with Level_*, aux and mage-tab
	if level == 15:
		#for pat in dir_include_patterns:
		auxSplit = targethost.split('aux')
		levelSplit = targethost.split('Level_')
		mageSplit = targethost.split('mage-tab')
		#matched = re.match('.*Level*', targethost) or re.match('.*aux*', targethost) or re.match('.*mage-tab*', targethost)
		if len(auxSplit) > 1 or len(levelSplit) > 1 or len(mageSplit) > 1:
			print 'Processing matched - Make dir: %s level %s' % (targethost, str(level))
			make_path(localpath)    	    	  		
			queue.put(DownloadDir(targethost, localpath))
	else:
		print 'Processing  dir: %s level %s' % (targethost, str(level))
		make_path(localpath)	
		queue.put(DownloadDir(targethost, localpath))
	return

def main(argv):
    global numberOfDirThreads, numberOfFileThreads, errorLogFile	
    if len(sys.argv) != 4:
        print 'Illegal arguments. Example of proper usage: tcga_download.py <config file> <tumor type> <errorLog>'
        sys.exit(-1)

    config = ConfigParser.RawConfigParser()
    config.read(sys.argv[1])

    subFilePath = sys.argv[2]
    errorLogFile = open(sys.argv[3], 'w')

    host = config.get("Connection", "host")
    host.rstrip('/')

    username = config.get("Connection", "username")
    password = config.get("Connection", "password")

    localPath = config.get("Local", "localPath")
    localPath.rstrip('/')

    numberOfDirThreads = config.getint("Local", "numberOfDirThreads")
    numberOfFileThreads = config.getint("Local", "numberOfFileThreads")

    print 'Mirroring data on host %s for tumor type %s start [%s]' % (host, subFilePath, time.strftime("%c"))
    print ('  from %s' % host)
    print ('    to %s' % localPath)

    print 'Starting Dir Thread Pool with %s dir threads' % (numberOfDirThreads)
    for i in range(numberOfDirThreads):
        ts = DownloadDirThread(username=username, password=password)
        ts.setName("Dir" + str(i))
        ts.setDaemon(True)
        ts.start()
 
    print 'Starting File Thread Pool with %s threads' % (numberOfFileThreads)
    for i in range(numberOfFileThreads):
        fs = DownloadFileThread(username=username, password=password)
        fs.setName("File" + str(i))
        fs.setDaemon(True)
        fs.start()
        
    handle_subdir(host + '/' + subFilePath, localPath + '/' + subFilePath)
    #print 'Mirroring data for tumor type %s completed [%s]' % (subFilePath, time.strftime("%c"))

    queue.join()
    fileQueue.join()

if __name__ == "__main__":
    main(sys.argv[1:])

