import os
import sys
import re
import shutil
import stat

class IsFileError(IOError):
    def __init__ (self, boguspath):
        self.fault = boguspath
    def __str__ (self):
        return "IsFileError implies that this \
is a file when I wanted a directory: " + self.fault
    def __repr__(self):
        return "IsFileError(" + repr(self.fault) +")"

class DestinationClobberError(IOError):
    def __init__ (self, boguspath):
        self.fault = boguspath
    def __str__ (self):
        return "DestinationClobberError implies that this is probably not\
 something you wanted me to delete: " + self.fault
    def __repr__(self):
        return "DestinationClobberError(" + repr(self.fault) +")"


level_regex = re.compile(r"(?P<key>.*\.Level_\d+\.\d+\.)"
                        r"(?P<revision>\d+)\.(?P<build>\d+)$"
                        )
magetab_regex = re.compile(r"(?P<key>.*\.mage-tab\.\d+\.)"
                        r"(?P<revision>\d+)\.(?P<build>\d+)$"
                        )

def max_shove(dic, key, value):
    """Update a key in the dictionary if the previous value was lower.
    
    dic -- dictionary to insert into
    key -- key to query
    value -- proposed new value
    
    The dictionary is queried for the key. If the key is not in the
    dictionary at all, it is inserted with the specified value.
    If it is in the dictionary, its corresponding value is compared
    to the provided value. If the new value is higher, the value in
    the dictionary is updated and the function returns True (same
    as if there is no key); if it is not, the dictionary is unchanged
    and the function returns False."""
    if (key not in dic) or value > dic[key]:
        dic[key] = value
        return True
    return False

def slay_the_weaker(dic, key, value):
    """Delete the less-up-to-date directory by the same name, if any.
    
    Parameters:
    dic -- Dictionary of absolute paths minus version numbers to
            version numbers
    key -- String representing the absolute path minus version number
            to check against
    value -- Version number of candidate for new maximum or deletion
    
    If the key does not exist in the table, it is created and this
    function returns False. If the key does exist, values are compared.
    The directory corresponding to the lower value is recursively
    deleted, and the higher value is inserted into the table. If the
    values are equal, the function simply returns, although this case
    shouldn't happen in the context of how the function was initially
    written."""
    
    if not(key in dic):
        dic[key] = value
        return False
    low = None
    #print "TWO DIRECTORIES ENTER, ONE DIRECTORY LEAVES"
    #print "==>    ", key, ">>>", value, dic[key]
    if value > dic[key]:
        low = key + str(dic[key][0]) + "." + str(dic[key][1])
        dic[key] = value
    elif value == dic[key]:
        print "Warning: slay_the_weaker found a clone:", key, ".", value
        return False;
    else:
        low = key + str(value[0]) + "." + str(value[1])
    print "Sweeping out of date directory", low
    shutil.rmtree(low)
    return True

def dont_kill_the_repository(path, recursive = True):
    """Throw an exception if the specified path is not safe to write.
    
    Parameters:
    path -- Path to check for safety
    recursive -- Whether or not to check subdirectories for safety,
                 as opposed to blindly assuming it
    Throws:
        A DestinationClobberError if it looks like the specified directory
        is anything other than a snapshot directory.
    A snapshot directory should contain only symlinks and subdirectories
    that are themselves valid snapshot directories. Anything else
    might be real data and is not safe to blindly overwrite with this
    specific tool, which should only be writing into a new directory or
    overwriting a previous snapshot directory because its write behavior
    is highly destructive."""
    
    if not os.path.isdir(path):
        raise DestinationClobberError(path)
    for name in os.listdir(path):
        checkpath = os.path.join(path, name)
        if not os.path.islink(checkpath):
            if os.path.isdir(checkpath):
                if recursive:
                    dont_kill_the_repository(checkpath)
            else:
                raise DestinationClobberError(checkpath)
    return True

def recursive_symlinkify(src, dest):
    """Create a recursive symlink copy in dest of src.
    
    Parameters:
    src -- Path to make a shallow symlink copy of
    dest -- Where to put it
    If dest is already populated, it is verified to be safe to overwrite
    using a recursive call to dont_kill_the_repository. It is then
    deleted and simply rewritten; no attempt is made to "update" such
    a directory. A directory is safe to overwrite if and only if it
    contains only symlinks and subdirectories that are, by the same
    rule, safe to overwrite."""
    if not os.path.isdir(src):
        raise IsFileError(src)
    if os.path.exists(dest):
        if os.path.isfile(dest):
            raise DestinationClobberError(destpath)
        if os.path.isdir(dest):
            dont_kill_the_repository(dest)
            shutil.rmtree(dest)
        else:
            raise IOError
    os.mkdir(dest)
    for name in os.listdir(src):
        srcpath = os.path.join(src, name)
        destpath = os.path.join(dest, name)
        if os.path.isdir(srcpath):
            recursive_symlinkify(srcpath, destpath)
        elif os.path.isfile(srcpath):
            os.symlink(srcpath, destpath)

def sweep(path):
    """Delete everything but the most up-to-date revisions."""
    dont_kill_the_repository(path, False)
    matchtable = {}
    for entry in os.listdir(path):
        fullpath = os.path.join(path, entry)
        if not os.path.isdir(fullpath):
            continue
        match = level_regex.match(entry)
        if match == None:
            match = magetab_regex.match(entry)
        if match == None:
            continue
        #now it's one of the types
        myrevision = (int(match.group("revision")),
                      int(match.group("build")))
        slay_the_weaker(matchtable, 
                        os.path.join(path, match.group("key")),
                        myrevision)

def recursive_latest(src, dest):
    """Create a snapshot of the latest data in the directory.
    
    Recurses on all subdirectories. Latest is defined as, for each
    major version of a MAGE_TAB directory, the highest-numbered minor
    version; for a leve_n directory, for each major.minor version,
    the highest-numbered revision. (Assumes major.minor.revis.build)"""
    
    #destination directory-making postponed: only happens if
    #there's anything in this entire tree that we care about
    
    matchtable = {}
    
    for entry in os.listdir(src):
        srcpath = os.path.join(src, entry)
        destpath = os.path.join(dest, entry)
        if os.path.isfile(srcpath):
            continue
        print "at", src, "processing", entry
        match = level_regex.match(entry)
        if match == None:
            match = magetab_regex.match(entry)
        #now if it's none, it's neither type
        if match != None:
            myrevision = (int(match.group("revision")), 
                          int(match.group("build")))
            max_shove(matchtable, match.group("key"), myrevision)
        else:
            recursive_latest(srcpath, destpath)
            
    #matchtable now contains the highest edition for each
    print matchtable
    if len(matchtable):
        if not os.path.exists(dest):
            os.makedirs(dest)
        elif os.path.isfile(dest):
            raise IsFileError(dest)
    for primary, edition in matchtable.iteritems():
        name = primary + str(edition[0]) + "." + str(edition[1])
        srcpath = os.path.join(src, name)
        destpath = os.path.join(dest, name)
        recursive_symlinkify(os.path.abspath(srcpath), destpath)
    
    if len(matchtable):
        sweep(dest)
    
    return len(matchtable)

def is_usage(arg):
    """Check if someone probably wanted the usage string."""
    return arg in ["-h", "-?", "-H", "--help", "--usage"]

def recursive_chmod(path, value):
    """Recursively sets directory permissions."""
    for root, dirs, files in os.walk(path):
        for directory in dirs:
            os.chmod(os.path.join(root, directory), value)

if __name__ == "__main__":
    if len(sys.argv) < 3 or is_usage(sys.argv[1]):
        print """Creates a snapshot of the latest data in a dcc tree.
        Data is identified as a directory ending in Level_[n].[version]
        or mage-tab.[version], where a Level directory is expected
        to have a four-part version and a mage-tab expected to have
        a three-part version.
        
        The "latest" data is the highest-revision Level data for each
        uniquely-named level_* root and each major.minor version, and
        the highest-minor-version mage-tab data for each uniquely-named
        mage-tab root.
        
        The snapshot is a clone of the directory structure
        stripped down to directories that contain data,
        and data directories themselves only have the latest, which
        is a symlink to the real data directory.
        
        USAGE: dcc-snapshot-script.py src dest
            src: the dcc directory to take a snapshot of
            dest: the directory to put the snapshot tree in
            
            Caution: no verification is performed that src and dest
            are sane or distinct. However, the order in which operations
            are performed guarantees that overlapping src and dest
            will not result in infinite recursion. No other guarantees
            are made in this case."""
        exit(0)
    recursive_latest(sys.argv[1], sys.argv[2])
    recursive_chmod(sys.argv[2],
        stat.S_IRWXU | stat.S_IRGRP | stat.S_IXGRP | stat.S_IROTH 
        | stat.S_IXOTH)
