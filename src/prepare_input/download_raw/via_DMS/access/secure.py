from os import environ
class Config:
    # Set the database credentials in the environment.
    db_user = environ.get("DATABASE_USERNAME")
    db_password = environ.get("DATABASE_PASSWORD")
    db_server = environ.get("DATABASE_SERVER")
    db_name = environ.get("DATABASE_NAME")
