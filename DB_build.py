if __name__ == '__main__':

    from DB_schema import Base
    from sqlalchemy import create_engine
    from sqlalchemy.orm import sessionmaker
    from Parsers.Parser import parse_all
    from DB_Statistics.DB_Statistics import get_db_stats_interactors, get_db_stats_interactions

    #('sqlite:///C:\\Users\\olgas\\Desktop\\PaIntDB.db')
    #'sqlite:////Users/olga/Desktop/PaIntDB.db'
    engine = create_engine('sqlite:////Users/olga/Desktop/PaIntDB.db')
    # if db is already created, next line can be commented out (tables do not have to be created again
    # to connect to the db)
    #Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    # comment out whatever processes are unnecessary
    #parse_all(session)
    # get_db_stats_interactors('PAO1',session)
    #get_db_stats_interactions('PAO1',session)
    # get_db_stats_interactors('PA14', session)
    #get_db_stats_interactions('PA14', session)
