#!/usr/bin/python

# NIH TCGA downloader

import httplib2
import os, sys, time, errno, re, csv

public_host = 'http://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/anonymous/tumor/'
secure_host = 'https://tcga-data.nci.nih.gov/tcgafiles/ftp_auth/distro_ftpusers/tcga4yeo/tumor/'

host = ''
local_base = ''#'/titan/cancerregulome3/TCGA/repositories/dcc-mattiannala'

user = ''#'ShmulevI'
password = ''# 'CancerReg2***'

file_index = {}

dir_skip_patterns = [
	'^minbio',
	'^bio$',
	'^tracerel$',
	'tissue_images',
	#'.*Level_2.*',
	#'.*Level_3.*',
	'.*\.aux\..*',
	'archive',
	#'hg-cgh-415k_g4124a',
	#'mskcc.org',
	'^gsc$',
	'humanhap',
	
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
	'.*\.tsv$',
	'.*\.gif$',
	'.*_genes\.txt$',
	'.*_tags\.txt$',
	'.*_frequency\.txt$',
	'nohup.out',
]



def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError, exc:
        if exc.errno == errno.EEXIST: pass
        else: raise



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
		m = re.match('(.+)\.Level_1\.(\d+)\.(\d+)\.\d+$', dir)
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



def file_size(conn, url):
	resp, content = conn.request(url, 'HEAD', 
		headers={'Accept-Encoding': 'plain'})
	return int(resp['content-length'])



def handle_subdir(conn, dir):
	time.sleep(1)
	print 'Procesing dir: %s' % (dir)

	resp, content = conn.request(host + dir, 'GET')
	if resp['status'] != '200':
		print resp
		print 'HTTP request failed for ' + dir 
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
		mkdir_p(local_base + dir)
				
		ftp_size = file_size(conn, host + dir + '/' + ffile)

		ffile_path = local_base + dir + '/' + ffile

		local_size = 0
		if not os.path.exists(ffile_path):
			if os.path.exists(ffile_path + '.bz2'):
				local_size = ftp_size
			elif os.path.exists(ffile_path + '.gz'):
				local_size = ftp_size
		else:
			local_size = os.path.getsize(ffile_path)
		
		if local_size != 0 and local_size == ftp_size:
			print '= %s/%s' % (dir, ffile)
			sys.stdout.flush()
			continue
		
		print '+ %s/%s' % (dir, ffile)
		sys.stdout.flush()
		
		if os.system('wget -O %s%s/%s -q --user=%s --password=%s %s%s/%s' % 
			(local_base, dir, ffile, user, password, host, dir, ffile)) != 0:
			sys.exit('Download failed.')

	for fsub in filter_subdirs(subdirs):
		handle_subdir(conn, dir + '/' + fsub )
	
	return




def read_file_index(index_path):
	file = csv.reader(open(index_path), delimiter='\t', quoting=csv.QUOTE_NONE)
	index = {}
	for row in file:
		filename = row[0]
		filename = filename[2:]
		filesize = int(row[1])
		index[filename] = filesize
	print index
	return index




def save_file_index(index, index_path):
	file = open(index_path, 'w')
	for (key, val) in index.iteritems():
		file.write('%s\t%d\n' % (key, val))
	file.close()



def main(argv):
	global host
	global file_index
	global user
	global password
	global local_base
	
	if len(sys.argv) != 5:
		print 'Bad amount of arguments. Example of proper usage:'
		print './tcga_download.py <tumor type> <local_base> <user> <password>'
		sys.exit(-1)
	
	timenow = time.strftime("%c")
	print 'Mirroring DCC begins: %s' % (timenow)

	tumor_type = sys.argv[1].lower()
	local_base = sys.argv[2]
	user = sys.argv[3]
	password = sys.argv[4]	
	
	#file_index = read_file_index(
	#	local_base + '/' + tumor_type + '/filelist.txt')
	
	print 'Mirroring public data for tumor type %s:' % tumor_type.upper()
	print (public_host + tumor_type)
	
	conn = httplib2.Http()
	host = public_host
	handle_subdir(conn, tumor_type)

	print
	print 'Mirroring restricted data for tumor type %s:' % tumor_type.upper()
	print (secure_host + tumor_type)

	conn = httplib2.Http()
	host = secure_host
	conn.add_credentials(user, password)
	handle_subdir(conn, tumor_type)
	
	#save_file_index(file_index,
	#	local_base + '/' + tumor_type + '/filelist.txt')
	timenow = time.strftime("%c")
	print 'Mirroring DCC completed: %s' % (timenow)
	return

if __name__ == "__main__":
	main(sys.argv[1:])

