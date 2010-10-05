import sys
import os
from time import gmtime, strftime

#import re
import ConfigParser

#mostRecentMsPattern = r'''(?P<base>.+)-(?P<version>[\w.]+)-(?P<release>[\w.]+)\.(?P<arch>\w+)\.rpm'''
#mostRecentMsPattern = r'''(?P<kvalue>.+).(?P<mvalue>[\w.]+).0'''

#def appendMostRecentMs(conn, url, listOfFiles):
#    resp, content = conn.request(url, 'GET')
#    if resp['status'] != '200':
#        print resp
#        print 'HTTP request failed for ' + url
#        return
#
#    for line in content.splitlines():
#        m = re.match('^\s*<a href=.+?>(.+)</a>.*', line)
#        if not m: continue
#        print "line=" + line
#        listOfFiles.append(m.group(1))
#
#    for filename in listOfFiles:
#        m_value = get_m_value(filename)
#        print filename + ":" + m_value

#def get_m_value(filename):
#    patternRe = re.compile(mostRecentMsPattern, re.VERBOSE)
#
#    m = patternRe.search(filename)
#    if m:
#        return m.group('mvalue')
#
#    return None

#def downloadMostRecentMs(ci, targetPath):
#    targetPath.lstrip('/')
#
#    listOfFiles = []
#    appendMostRecentMs(ci.conn, ci.host + "/" + targetPath, listOfFiles)
#
#    for filename in listOfFiles:
#        print "filename=" + filename

def scan_all_files(path, list_of_files):
    path = remove_trailing_slash(path)

    if os.path.isfile(path):
        list_of_files.append(path)
    else:
        items = os.listdir(path)
        for item in items:
            scan_all_files(path + "/" + item, list_of_files)

def filter_latest(list_of_files):
    print "filter_latest(" + str(len(list_of_files)) + ")"

    # TODO : check paths against rules
    latest_files = []
    for file_path in list_of_files:
        if file_path.count("ov/cgcc/hms.harvard.edu/hg-cgh-244a/cna") == 1:
            latest_files.append(file_path)

    return latest_files

def generate_links(local_path, snapshot_path, latest_files):
    print "generate_links(" + local_path + "," + snapshot_path + "," + str(len(latest_files)) + ")"

    for file_path in latest_files:
        new_file_path = file_path.replace(local_path, snapshot_path)
        new_file_dir = os.path.dirname(new_file_path)
        if not os.path.exists(new_file_dir):
            os.makedirs(new_file_dir)

        os.symlink(file_path,new_file_path)

def remove_trailing_slash(target_string):
    if target_string is None:
        return target_string

    stripped = target_string.rstrip()
    if (stripped.endswith("/")):
        return stripped[:-1]

    return stripped

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print 'Illegal arguments. Example of proper usage: tcga_snapshot.py <config file>'
        sys.exit(-1)

    config = ConfigParser.RawConfigParser()
    config.read(sys.argv[1])

    local_path = remove_trailing_slash(config.get("Local", "localPath"))
    snapshot_path = remove_trailing_slash(config.get("Local", "snapshotPath"))

    list_of_files = []

    scan_all_files(local_path, list_of_files)

    latest_files = filter_latest(list_of_files)

    snapshot_dir = snapshot_path + "/" + strftime("%Y_%m_%d_%H%M%S", gmtime())
    print "creating snapshot directory:" + snapshot_dir
    os.mkdir(snapshot_dir)

    generate_links(local_path, snapshot_dir, latest_files)

#    downloadMostRecentMs(ci, "ov/cgcc/hms.harvard.edu/hg-cgh-244a/cna")
#    downloadMostRecentMs(ci, "ov/cgcc/hms.harvard.edu/hg-cgh-415k_g4124a/cna")
#    downloadMostRecentMs(ci, "ov/cgcc/broad.mit.edu/ht_hg-u133a/transcriptome")
#    downloadMostRecentMs(ci, "ov/cgcc/jhu-usc.edu/humanmethylation27/methylation")
#    downloadMostRecentMs(ci, "ov/cgcc/mskcc.org/cgh-1x1m_g4447a/cna")
#    downloadMostRecentMs(ci, "ov/cgcc/unc.edu/agilentg4502a_07_2/transcriptome")

#    target_dirs.append("ov/cgcc/hms.harvard.edu/illuminaga_mrna_dge/transcriptome")
#    scanPaths(config, mostRecentM_onlyK1("ov/cgcc/unc.edu/agilentg4502a_07_3/transcriptome"))
#    scanPaths(config, mostRecentM_onlyK1("ov/cgcc/unc.edu/h-mirna_8x15kv2/mirna"))

