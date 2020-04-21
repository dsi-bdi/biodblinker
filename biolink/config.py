from os.path import join, isfile, isdir
import gzip
import bz2
import configparser
import biolink
import pkg_resources
import codecs
import sys
import os

def file_open(filepath, *args, ** kwargs):
    """ Open a file (plain or compressed)

    Parameters
    ----------
    filepath : str
        file path to a data map file

    Returns
    -------
    file_descriptor
        file descriptor object
    """
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode="rt", *args, **kwargs)

    elif filepath.endswith('.bz2'):
        return bz2.open(filepath, mode="rt", *args, **kwargs)
    else:
        return open(filepath, mode="rt", *args, **kwargs)


def load_data_map(filepath):
    """ Load  a data mapping file

    Parameters
    ----------
    filepath : str
        file path to a data map file

    Returns
    -------
    dict
        mapping dictionary
    """
    data_dict = dict()
    fqfp = join(get_data_directory(), filepath)
    if isfile(fqfp):
        with file_open(fqfp) as data_fd:
            for line in data_fd:
                key, val = line.strip().split("\t")
                if val == '-':
                    data_dict[key] = []
                else:
                    vals = val.split(';')
                    data_dict[key] = vals

    else:
        raise FileNotFoundError(f"The file ({fqfp}) does not exist.")

    return data_dict


def get_database_mapping_sources(database_name):
    """ Get database mapping files by database name

    Parameters
    ----------
    database_name : str
        database name

    Returns
    -------
    dict
        file sources key names and paths
    """
    config = configparser.ConfigParser()

    with pkg_resources.resource_stream(biolink.__name__, 'sources.ini') as stream:
        config.read_file(codecs.getreader("utf-8")(stream))
        if database_name in set(config.sections()):
            mapping_dict = dict(config[database_name])
            return {map_name: map_short_path for map_name, map_short_path in mapping_dict.items()}
        else:
            raise ValueError(f"Database name ({database_name}) could not be found in the source.ini file.")


_data_dir = None


def get_data_directory():
    global _data_dir
    if _data_dir is not None:
        return _data_dir

    if 'BIOLINK' in os.environ:
        data_dir = os.environ['BIOLINK']
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        elif not os.path.isdir(data_dir):
            raise ValueError(f'BIOLINK environment variable is a file not a directory')
    elif sys.platform == 'win32' and 'APPDATA' in os.environ:
        data_dir = os.environ['APPDATA']
        data_dir = os.path.join(data_dir, 'biolink_data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
    else:
        data_dir = os.path.expanduser('~')
        data_dir = os.path.join(data_dir, 'biolink_data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)

    _data_dir = data_dir
    return data_dir

