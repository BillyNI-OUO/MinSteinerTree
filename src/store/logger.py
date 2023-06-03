import logging
from ..constants import Constants
class CustomFormatter(logging.Formatter):

    grey = "\x1b[38;20m"
    yellow = "\x1b[33;20m"
    red = "\x1b[31;20m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(levelname)s - %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)
    
class storeLogger:
    logger:logging.Logger
    
    def __init__(self, filename:str=Constants.LOG_FILE_NAME, level=Constants.LOG_LEVEL) -> None:
        
        logging.basicConfig(level=level, format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p', force=True)
        logging.getLogger('matplotlib.font_manager').disabled = True
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        self.info(f'Setting up Logger')

        
    def info(self, msg:str)->None:
        msg = f'{Constants.GREY_STYLE} [INFO]: {msg}{Constants.RESET_STYLE}'
        logging.info(msg=msg)

    def debug(self, msg:str)->None:
        msg = f'{Constants.GREY_STYLE} [DEBUG]: {msg}{Constants.RESET_STYLE}'
        logging.debug(msg=msg)

    def error(self, msg:str)->None:
        msg = f'{Constants.RED_STYLE} [ERROR]: {msg}{Constants.RESET_STYLE}'
        logging.error(msg=msg)




Logger = storeLogger()
