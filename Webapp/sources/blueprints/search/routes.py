from flask import Response, jsonify, render_template, request
from flask_login import login_required
from sources import models, services
from sources.auxiliary import get_workbooks, get_workgroups, smiles_to_inchi
from sources.extensions import db

from . import search_bp


@search_bp.route("/updated_workgroup_dropdown", methods=["POST"])
@login_required
def updated_workgroup_dropdown() -> Response:
    """when the workgroup dropdown is updated, gets all workbooks which belong to the selected workgroup"""
    # institution = request.form['institution']
    workgroup = request.form["workgroup"]
    workbooks = get_workbooks(workgroup)
    workbooks.insert(0, "All")
    return jsonify({"workbooks": workbooks})


# @search_bp.route("/text_search_handler", methods=["POST"])
# @login_required
# def text_search_handler() -> Response:
#     """Gets search options for text.
#     **not yet implemented**
#     """
#     workgroup = request.form["workgroup"]
#     workbook = request.form["workbook"]


@search_bp.route("/structure_search_handler", methods=["POST"])
@login_required
def structure_search_handler() -> Response:
    """Receives ajax request from search.js"""
    search = SearchHandler()
    if not search.reactions:
        # if no reactions exit search and post exception to frontend
        return search.exit()
    # get the search type - only exact_structure currently.
    search_type = request.form["searchType"]
    if search_type == "exact_structure":
        return search.exact_structure_search()


class SearchHandler:
    def __init__(self):
        """On initialisation creates the list 'reactions'"""
        print("making search handler")
        self.workgroup = [request.form["workgroup"]]
        self.workbook = [request.form["workbook"]]
        self.reactions = []
        self.matches = []
        self.search_results = ""
        self.schemes = []
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
        target_inchi = smiles_to_inchi(smiles)
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
                    "schemes": self.schemes,
                }
            )
        else:
            return jsonify({"status": "fail", "message": "No results found"})

    def exact_structure_match_loop(self, reaction, target_inchi) -> None:
        """Exact structure search looking for matches in reactants, reagents, and products"""
        if reaction.reaction_smiles:
            reactants1, products1 = reaction.reaction_smiles.split(">>")
            reactants1, products1 = reactants1.split("."), products1.split(".")
            reactants2, reagents, products2 = (
                reaction.reactants,
                reaction.reagents,
                reaction.products,
            )
            reactants, products = list(set(reactants1 + reactants2)), list(
                set(products1 + products2)
            )
            components = [reactants, reagents, products]
            for component_type in components:
                for component_smiles in component_type:
                    if component_smiles:
                        inchi = smiles_to_inchi(component_smiles)
                        if inchi == target_inchi:
                            self.matches.append(reaction)
                            return

    def render_results_template(self) -> None:
        """Makes the results list"""
        reactions = services.reaction.to_dict(self.matches, "AZ")
        self.schemes = services.reaction.make_scheme_list(self.matches, "small")
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
