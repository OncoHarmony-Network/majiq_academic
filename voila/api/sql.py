import os
from abc import abstractmethod, ABC

from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker

Session = sessionmaker()


class SQLType(ABC):
    @abstractmethod
    def __init__(self, session):
        self.session = session

    def __bool__(self):
        return self.exists

    def __iter__(self):
        return self.get.__iter__()

    @abstractmethod
    def add(self, **kwargs):
        pass

    @property
    @abstractmethod
    def get(self):
        pass

    @property
    @abstractmethod
    def exists(self):
        pass


class SQL:
    def __init__(self, filename, model, delete=False, ):
        if delete is True:
            try:
                os.remove(filename)
            except FileNotFoundError:
                pass

        engine = create_engine('sqlite:///{0}'.format(filename))
        event.listen(engine, 'connect', self._fk_pragma_on_connect)
        model.Base.metadata.create_all(engine)
        Session.configure(bind=engine)

        self.session = Session()
        self.filename = filename

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @staticmethod
    def _fk_pragma_on_connect(dbapi_con, con_record):
        dbapi_con.execute('pragma foreign_keys=ON')

    def commit(self):
        self.session.commit()

    def close(self):
        self.commit()
        self.session.close_all()
