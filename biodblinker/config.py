from os.path import join, isfile, isdir, exists
import gzip
import bz2
import configparser
import biodblinker
import pkg_resources
import codecs
import sys
import os
from zipfile import ZipFile
from os import remove
from biodblinker.fileio import download_file, read_remote_checksum

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

    with pkg_resources.resource_stream(biodblinker.__name__, 'sources.ini') as stream:
        config.read_file(codecs.getreader("utf-8")(stream))
        if database_name in set(config.sections()):
            mapping_dict = dict(config[database_name])
            return {map_name: map_short_path for map_name, map_short_path in mapping_dict.items()}
        else:
            raise ValueError(f"Database name ({database_name}) could not be found in the source.ini file.")


def get_all_mappings_sources():
    """ Get all database mapping files

    Parameters
    ----------

    Returns
    -------
    dict
        file sources key names and paths
    """
    config = configparser.ConfigParser()
    mapping_sources = {}
    with pkg_resources.resource_stream(biodblinker.__name__, 'sources.ini') as stream:
        config.read_file(codecs.getreader("utf-8")(stream))
        for section in set(config.sections()):
            if section == 'default':
                continue
            mapping_dict = dict(config[section])
            mapping_sources[section] = {map_name: map_short_path for map_name, map_short_path in mapping_dict.items()}
    return mapping_sources

_data_dir = None


def get_data_directory():
    global _data_dir
    if _data_dir is not None:
        return _data_dir

    if 'biodblinker' in os.environ:
        data_dir = os.environ['biodblinker']
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
        elif not os.path.isdir(data_dir):
            raise ValueError(f'biodblinker environment variable is a file not a directory')
    elif sys.platform == 'win32' and 'APPDATA' in os.environ:
        data_dir = os.environ['APPDATA']
        data_dir = os.path.join(data_dir, 'biodblinker_data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
    else:
        data_dir = os.path.expanduser('~')
        data_dir = os.path.join(data_dir, 'biodblinker_data')
        if not os.path.exists(data_dir):
            os.mkdir(data_dir)

    _data_dir = data_dir
    return data_dir


def verify_mappings():
    all_exist = True
    biodblinker_data = get_data_directory()
    mapping_sources = get_all_mappings_sources()
    for _, mappings in mapping_sources.items():
        for map_name, short_path in mappings.items():
            map_path = join(biodblinker_data, short_path)
            valid_map = exists(map_path) & isfile(map_path)
            if not valid_map:
                all_exist = False
                break

    return all_exist


def download_mappings():
    biodblinker_data = get_data_directory()
    checksum_uri = f'https://github.com/dsi-bdi/biodblinker_mappings/releases/download/v{biodblinker.__data_version__}/biodblinker_mappings.zip.md5'
    checksum = read_remote_checksum(checksum_uri)
    mappings_uri = f'https://github.com/dsi-bdi/biodblinker_mappings/releases/download/v{biodblinker.__data_version__}/biodblinker_mappings.zip'

    mappings_fp = join(biodblinker_data, "./biodblinker_mappings.zip")
    download_file(mappings_uri, mappings_fp, checksum)

    with ZipFile(mappings_fp, 'r') as map_zip:
        map_zip.extractall(path=biodblinker_data)

    os.remove(mappings_fp)


_config_verified = False


def verify_biodblinker():
    global _config_verified
    if _config_verified:
        return
    if not verify_mappings():
        print('biodblinker mappings do not exist, downloading mappings size 128MB')
        download_mappings()
        _config_verified = True
