import logging, functools, time, datetime
logger = logging.getLogger('decorator-log')
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)
# disable logging on the standard error stream
logger.disabled = False  

class LogDecorator(object):
    def __init__(self):
        self.logger = logging.getLogger('decorator-log')
        self.log_debug_file     = 'log_debug.dat'
        self.log_exception_file = 'log_exception.txt'
        self.debug = False

    def __call__(self, fn):
        @functools.wraps(fn)
        def decorated(*args, **kwargs):
            try:
                # .{fn.__name__}
                msj = f'>> {fn.__module__}.{fn.__qualname__}(): >> called with >>{args} - {kwargs}'
                self.logger.debug(msj)
                self.append_log_to_file(self.log_debug_file, msj )
                result = fn(*args, **kwargs)
                self.logger.debug(result)
                return result
            except Exception as ex:
                msj = f"Exception raised in >> {fn.__module__}.{fn.__qualname__}(): >> exception: {str(ex)}\n"
                self.logger.debug(msj)
                self.append_log_to_file(self.log_exception_file, msj )
                raise ex
            return result
        return decorated

    def append_log_to_file(self, file, message: str) -> None:
        log_file = open(file, 'a')
        log_file.write(f'>> {datetime.datetime.now()} >> {message}\n')
        log_file.close()

    def reset_log_file(self):
        log_file = open(file, 'w+')
        log_file.write('')
        log_file.close()
        