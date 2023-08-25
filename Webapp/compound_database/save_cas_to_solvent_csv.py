import Compound_database_extraction as CDE
from config import basedir
import os
from pony.orm import select, db_session


def link_solvent_and_compound():
	# connect to db
	db_config = {
		'provider': 'sqlite',
		'filename': os.environ.get('DATABASE_URL') or os.path.join(basedir, 'test_database.sqlite'),
		'create_db': False
	}
	db = CDE.open_database(db_config)
	# iterate through solvents
	with db_session:
		solvent_ls = select(x for x in db.Solvent)[:]
		for sol in solvent_ls:
			corresponding_compound_db_entry = select(x for x in db.Compound if x.name.lower() == sol.name.lower()).first()
			sol.set(**{'compound': corresponding_compound_db_entry})


if __name__ == '__main__':
	link_solvent_and_compound()