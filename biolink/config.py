from os.path import join, isfile, isdir
import gzip
import bz2
import configparser
import biolink


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
    if isfile(filepath):
        with file_open(filepath) as data_fd:
            for line in data_fd:
                key, val = line.strip().split("\t")
                if val == '-':
                    data_dict[key] = []
                else:
                    vals = val.split(';')
                    data_dict[key] = vals

    else:
        raise FileNotFoundError(f"The file ({filepath}) does not exist.")

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
    lib_root = biolink.__LIB_ROOT__DIR__
    config = configparser.ConfigParser()
    config.read(join(lib_root, "sources.ini"))
    if database_name in set(config.sections()):
        mapping_dict = dict(config[database_name])
        return {map_name: join(lib_root, map_short_path) for map_name, map_short_path in mapping_dict.items()}
    else:
        raise ValueError(f"Database name ({database_name}) could not be found in the source.ini file.")
