if __name__ == '__main__':

    from Schema1 import Base
    from sqlalchemy import create_engine, or_
    from sqlalchemy.orm import sessionmaker

    engine = create_engine('sqlite:///:memory:', echo=True)
    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()


    def find_type(psi_code):
        return {
            '1110': 'predicted interaction',
            '2232': 'molecular association',
            '0914': 'association',
            '0915': 'physical association',
            '0407': 'direct interaction',
        }[psi_code]

    from Data.file_desc import parse

    parse()

