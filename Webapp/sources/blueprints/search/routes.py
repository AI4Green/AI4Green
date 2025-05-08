from json.decoder import JSONDecodeError

from flask import Response, json, jsonify, render_template, request
from flask_login import login_required
from sources import models, services
from sources.auxiliary import get_workbooks, get_workgroups
from sources.extensions import db

from . import search_bp

# @search_bp.route("/text_search_handler", methods=["POST"])
# @login_required
# def text_search_handler() -> Response:
#     """Gets search options for text.
#     **not yet implemented**
#     """


@search_bp.route("/structure_search_handler", methods=["POST"])
@login_required
@search_bp.doc(security="sessionAuth")
def structure_search_handler() -> Response:
    """
    Returns search results for exact structure search.

    Returns:
        flask.Response: returns the search results as a json object
    """
    search = SearchHandler()
    if not search.reactions:
        # if no reactions exit search and post exception to frontend
        return search.exit()
    # get the search type - only exact_structure currently.
    search_type = request.form["searchType"]
    if search_type == "exact_structure":
        mol = request.form["mol"]
        if "SRU" in mol:  # polymer found
            return search.exact_polymer_structure_search()
        else:
            return search.exact_structure_search()


class SearchHandler:
    def __init__(self):
        """On initialisation creates the list 'reactions'"""
        self.workgroup = [request.form["workgroup"]]
        self.workbook = [request.form["workbook"]]
        self.reactions = []
        self.matches = []
        self.search_results = ""
        self.images = []
        self.get_search_workgroups()
        self.get_search_workbooks()
        self.get_reactions()

    def get_search_workgroups(self):
        """Gets the workgroups that user has chosen for searching"""
        # if all workgroups - else if a chosen workgroup can use the user input as is
        if self.workgroup == ["All"]:
            self.workgroup = list(get_workgroups())

    def get_search_workbooks(self):
        """Gets the workbooks that user has chosen for searching"""
        # if all workbooks - else if a chosen workbook can use the user input as is
        if self.workbook == ["All"]:
            self.workbook = []
            self.workbook.extend(
                get_workbooks(workgroup) for workgroup in self.workgroup
            )
            # flatten nested lists
            self.workbook = [
                workbook for sublist in self.workbook for workbook in sublist
            ]

    def get_reactions(self):
        """Gets the reactions that user has chosen for searching"""
        # to get reactions - need to get list of workbooks
        for workbook in self.workbook:
            reactions = (
                db.session.query(models.Reaction)
                .filter(models.Reaction.status == "active")
                .join(models.WorkBook)
                .filter(models.WorkBook.name == workbook)
                .join(models.WorkGroup)
                .filter(models.WorkGroup.name.in_(self.workgroup))
                .all()
            )
            self.reactions.append(reactions)
        # flatten list
        self.reactions = [
            reaction for sublist in self.reactions for reaction in sublist
        ]

    def exact_structure_search(self) -> Response:
        """Performs an exact structure search on the reaction list"""
        # get search smiles and convert to inchi
        smiles = request.form["smiles"]
        target_inchi = services.all_compounds.smiles_to_inchi(smiles)
        if not target_inchi:
            return jsonify(
                {"status": "fail", "message": "Invalid structure, please try again!"}
            )
        # iterate through reactant, reagent, and product for each reaction and check for a match
        for reaction in self.reactions:
            self.exact_structure_match_loop(reaction, target_inchi)
        self.sort_matches()
        if self.matches:
            self.render_results_template()
            return jsonify(
                {
                    "status": "success",
                    "message": f"{len(self.matches)} results found",
                    "search_results": self.search_results,
                    "images": self.images,
                }
            )
        else:
            return jsonify({"status": "fail", "message": "No results found"})

    def exact_structure_match_loop(
        self, reaction: models.Reaction, target_inchi: str
    ) -> None:
        """Exact structure search looking for matches in reactants, reagents, and products"""
        if reaction.reaction_smiles:
            (
                reactants_from_reaction_smiles,
                products_from_reaction_smiles,
            ) = reaction.reaction_smiles.split(">>")
            (
                reactants_from_reaction_smiles,
                products_from_reaction_smiles,
            ) = reactants_from_reaction_smiles.split(
                "."
            ), products_from_reaction_smiles.split(
                "."
            )
            reactants_from_db, reagents, products_from_db = (
                reaction.reactants,
                reaction.reagents,
                reaction.products,
            )
            reactants, products = list(
                set(reactants_from_reaction_smiles + reactants_from_db)
            ), list(set(products_from_reaction_smiles + products_from_db))
            components = [reactants, reagents, products]
            for component_type in components:
                for component_smiles in component_type:
                    if component_smiles:
                        inchi = services.all_compounds.smiles_to_inchi(component_smiles)
                        if inchi == target_inchi:
                            self.matches.append(reaction)
                            return

    def exact_polymer_structure_search(self) -> Response:  # different func if polymer
        """Performs an exact structure search on the reaction list"""
        # get search smiles
        smiles = request.form["smiles"]
        smiles = smiles.split(" |")[0]
        if smiles.count("*") == 2:
            smiles = services.polymer_novel_compound.canonicalise(smiles)
        else:
            return jsonify(
                {
                    "status": "error message",
                    "message": "Searching for polymers with known end groups is not yet supported",
                }
            )  # TODO: deal with end groups and dummy atoms somehow

        # iterate through reactant, reagent, and product for each reaction and check for a match
        for reaction in self.reactions:
            self.exact_polymer_structure_match_loop(reaction, smiles)
        self.sort_matches()
        if self.matches:
            self.render_results_template()
            return jsonify(
                {
                    "status": "success",
                    "message": f"{len(self.matches)} results found",
                    "search_results": self.search_results,
                    "images": self.images,
                }
            )
        else:
            return jsonify({"status": "fail", "message": "No results found"})

    def exact_polymer_structure_match_loop(
        self, reaction: models.Reaction, target_smiles: str
    ) -> None:
        """Exact structure search looking for matches in reactants, reagents, and products"""
        if not reaction.reaction_type.value == "POLYMER":  # skip non-polymer reactions
            return

        if reaction.reaction_smiles:
            (
                reactants_from_reaction_smiles,
                products_from_reaction_smiles,
            ) = reaction.reaction_smiles.split(">>")
            (
                reactants_from_reaction_smiles,
                products_from_reaction_smiles,
            ) = reactants_from_reaction_smiles.split(
                "."
            ), products_from_reaction_smiles.split(
                "."
            )
            reactants_from_db, reagents, products_from_db = (
                reaction.reactants,
                reaction.reagents,
                reaction.products,
            )

            # load json strings and flatten lists
            for smiles_list in reactants_from_db, products_from_db:
                for i in range(len(smiles_list)):
                    try:
                        smiles_list[i] = json.loads(smiles_list[i])
                    except JSONDecodeError:
                        smiles_list[i] = [smiles_list[i]]

            reactants, products = list(
                set(reactants_from_reaction_smiles + sum(reactants_from_db, []))
            ), list(set(products_from_reaction_smiles + sum(products_from_db, [])))
            components = [reactants, reagents, products]
            for component_type in components:
                for component_smiles in component_type:
                    if component_smiles:
                        component_smiles = component_smiles.split(" |")[0]
                        component_smiles = (
                            services.polymer_novel_compound.clean_polymer_smiles(
                                component_smiles
                            )
                        )
                        if component_smiles == target_smiles:
                            self.matches.append(reaction)
                            return

    def render_results_template(self) -> None:
        """Makes the results list"""
        reactions = services.reaction.to_dict(self.matches)
        self.images = services.reaction.make_reaction_image_list(self.matches)
        self.search_results = render_template(
            "_saved_reactions.html", reactions=reactions, sort_crit="AZ"
        )

    def sort_matches(self) -> None:
        """Sorts matches"""
        self.matches.sort(key=lambda reaction: reaction.name)

    @staticmethod
    def exit() -> Response:
        """Exits the search if there are no reactions"""
        print("exit")
        return jsonify(
            {"status": "fail", "message": "No reactions found to search through"}
        )
