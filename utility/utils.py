import time
import argparse
from datetime import datetime, timezone
import tracemalloc

from functools import wraps
from inspect import getcallargs, signature
from collections import OrderedDict, Iterable
from itertools import *
import six


import logging
import json_log_formatter


json_handler = logging.FileHandler(filename='./log/{}.log'.format(datetime.utcnow().strftime("%Y_%m_%d-%I_%M_%S_%p")) )
json_handler.setFormatter(json_log_formatter.JSONFormatter())
# json_handler.setFormatter(json_log_formatter.VerboseJSONFormatter())
logger = logging.getLogger('metaPro')
logger.addHandler(json_handler)
logger.setLevel(logging.INFO)

logger.info('NMDC-PNNL Meta-Proteomics workflow run record.')


def stats(method):
    def checked(*args, **kw):
        ts = time.time()
        tracemalloc.start()
        result = method(*args, **kw)
        current, peak = tracemalloc.get_traced_memory()
        te = time.time()
        logger.info(method.__name__, extra={'runtime': time.strftime("%H:%M:%S", time.gmtime((te - ts))),
                                            'current_m/m(MB)': current / 10 ** 6,
                                            'peak_m/m(MB)'   : peak / 10 ** 6})
        tracemalloc.stop()
        return result
    return checked


def timeit(method):
    def timed(*args, **kw):
        ts = time.time()
        result = method(*args, **kw)
        te = time.time()
        print("{} took {}.".format(method.__name__,time.strftime("%H:%M:%S" , time.gmtime((te - ts))) ))
        logging.info(msg= "{} took {}.".format(method.__name__,time.strftime("%H:%M:%S" , time.gmtime((te - ts))) ) )
        return result
    return timed

def str2bool(v):
    '''
    chnages userInput to a yes/no
    :param v: string
    :return: bool
    '''
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError("Boolean value expected. ['yes', 'true', 't', 'y', '1', 'no', 'false', 'f', 'n', '0']" )

def current_local_datetime():
    '''

    :return: current local date and time
    '''
    return datetime.now()

def current_UTC_datetime():
    '''

    :return: current UTC date and time.
    '''
    return datetime.now(timezone.utc)

# This file defines a decorator '@log_to()' that logs every call to a
# function, along with the arguments that function was called with. It
# takes a logging function, which is any function that accepts a
# string and does something with it. A good choice is the debug
# function from the logging module. A second decorator '@logdebug' is
# provided that uses 'logging.debug' as the logger.


def flatten(l):
    """Flatten a list (or other iterable) recursively"""
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, six.string_types):
            for sub in flatten(el):
                yield sub
        else:
            yield el

def getargnames(func):
    """Return an iterator over all arg names, including nested arg names and varargs.
    Goes in the order of the functions argspec, with varargs and
    keyword args last if present."""
    (argnames, varargname, kwargname, _) = signature(func)
    return chain(flatten(argnames), filter(None, [varargname, kwargname]))

def getcallargs_ordered(func, *args, **kwargs):
    """Return an OrderedDict of all arguments to a function.
    Items are ordered by the function's argspec."""
    argdict = getcallargs(func, *args, **kwargs)
    return OrderedDict((name, argdict[name]) for name in getargnames(func))

def describe_call(func, *args, **kwargs):
    yield "Calling %s with args:" % func.__name__
    for argname, argvalue in getcallargs_ordered(func, *args, **kwargs).items():
        yield "\t%s = %s" % (argname, repr(argvalue))

def log_to(logger_func):
    """A decorator to log every call to function (function name and arg values).
    logger_func should be a function that accepts a string and logs it
    somewhere. The default is logging.debug.
    If logger_func is None, then the resulting decorator does nothing.
    This is much more efficient than providing a no-op logger
    function: @log_to(lambda x: None).
    """
    if logger_func is not None:
        def decorator(func):
            @wraps(func)
            def wrapper(*args, **kwargs):
                for line in describe_call(func, *args, **kwargs):
                    logger_func(line)
                return func(*args, **kwargs)
            return wrapper
    else:
        decorator = lambda x: x
    return decorator

logdebug = log_to(logging.debug)