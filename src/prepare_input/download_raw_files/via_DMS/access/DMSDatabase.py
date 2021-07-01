import pymssql
import logging
import sys
from utility.utils import logger

class DMSDatabase:
    ''' Database connection class'''
    def __init__(self, config):
        '''

        :param config:
        '''
        self.SERVER = config.db_server
        self.USER = config.db_user
        self.PASSWORD = config.db_password
        self.DATABASE_NAME = config.db_name
        self.conn = None

    def open_connection(self):
        '''Connection to DMS MS SQLserver.
        '''
        try:
            if self.conn is None:
                self.conn = pymssql.connect(server = self.SERVER,
                                            user = self.USER,
                                            password = self.PASSWORD,
                                            database= self.DATABASE_NAME)
                logger.info(msg="CONNECTION: {}".format(self.conn))
        except pymssql.MySQLError as conn_err:
            logger.error(msg="SQLserver connection FAILD\n{}".format(self.SERVER,conn_err) )
            # sys.exit()
        finally:
            logger.info("SQL connection to {}:{}:{} opened successfully!".format(self.SERVER,self.DATABASE_NAME, self.USER))

    def run_query(self, query):
        '''Execute SQL query.
        :param query: SQL query
        :return: database cursor, a pointer to a specific rows within a query result.
        '''
        try:
            self.open_connection()
            cursor = self.conn.cursor(as_dict=True)
            cursor.execute(query) # generator object.
            return cursor
        except Exception as e:
            logger.error(msg="CURSOR FAILED:".format(e))