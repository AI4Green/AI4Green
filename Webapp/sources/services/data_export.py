import json
from typing import Dict, List, Optional

import chython.files
from sources import models, services


class ReactionDataFile:
    def __init__(self, db_reaction: models.Reaction, filename: str):
        """
        Filename should not include extension.
        """
        self.db_reaction = db_reaction
        self.reaction_container = self.make_reaction_container()
        self.metadata = ReactionMetaData(db_reaction).get_dict()
        self.filename = filename
        if self.reaction_container:
            self.reaction_container.meta.update(self.metadata)

    def make_reaction_container(self):
        """
        Makes a reaction container by using the SMILES obtained from the db reaction object.
        Full reaction SMILES including reactants, agents (aka reagents and solvents), and products
        """
        reaction_smiles = self.make_reaction_smiles()
        # return None if we cannot make a reaction_container from smiles - likely to be due to no smiles present.
        try:
            return chython.files.smiles(reaction_smiles)
        except ValueError:
            return None

    def make_reaction_smiles(self) -> str:
        """
        Makes a reaction smiles in format reactant1.reactants2>reagents.solvents>products from a db reaction object

        Returns:
            the reaction smiles string.
        """
        reactant_smiles = remove_default_data(self.db_reaction.reactants)
        reagent_smiles = remove_default_data(self.db_reaction.reagents)
        solvent_smiles = services.all_compounds.get_smiles_list(
            self.db_reaction.solvent
        )
        product_smiles = remove_default_data(self.db_reaction.products)
        reaction_smiles = (
            ".".join(reactant_smiles)
            + ">"
            + ".".join(reagent_smiles)
            + ".".join(solvent_smiles)
            + ">"
            + ".".join(product_smiles)
        )
        return reaction_smiles

    def save(self):
        if self.reaction_container:
            self.filename += ".rdf"
            self.save_as_rdf()
        else:
            self.filename += ".json"
            self.save_as_json()

    def save_as_rdf(self):
        with chython.files.RDFWrite(self.filename) as f:  # context manager supported
            # for r in data:
            f.write(self.reaction_container)
        #  TESTING CODE @
        # opening file for testing purposes
        with chython.files.RDFRead(self.filename) as f:
            rdf_contents = next(f)

        non_string_types_inside_meta = [
            "reagents",
            "reactants",
            "solvents",
            "products",
            "standard_protocols_used",
            "file_attachment_names",
            "addenda",
            "solvent_sustainability",
            "sustainability_data",
        ]
        import ast

        for key in self.metadata.keys():
            # if the string is equal to none we reload
            if (
                rdf_contents.meta[key] == "None"
                or key in non_string_types_inside_meta
                and isinstance(rdf_contents.meta[key], str)
            ):
                print(key)
                rdf_contents.meta.update(
                    {key: ast.literal_eval(rdf_contents.meta[key])}
                )

            if self.reaction_container.meta[key] != rdf_contents.meta[key]:
                print(
                    key,
                    rdf_contents.meta[key],
                    self.reaction_container.meta[key],
                    "unequal",
                )

            if self.reaction_container.meta[key] == rdf_contents.meta[key]:
                print("equal")

        if (
            not self.reaction_container.meta["products"]
            == rdf_contents.meta["products"]
        ):
            print("prod")

        if (
            not self.reaction_container.meta["reagents"]
            == rdf_contents.meta["reagents"]
        ):
            print("reagents")

        if (
            not self.reaction_container.meta["reactants"]
            == rdf_contents.meta["reactants"]
        ):
            print("reactants")

        assert (
            rdf_contents == self.reaction_container
        ), "change in data during file read/write"

        print("we all good")

    def save_as_json(self):
        with open(self.filename, "w") as f:
            json.dump(self.metadata, f)

        with open(self.filename, "r") as f:
            meta_reload = json.load(f)

        assert meta_reload == self.metadata, "change in data during file read/write"


class ReactionMetaData:
    """
    Takes a reaction object from the database and makes a dictionary for metadata to describe the reaction
    For consistency, all missing or default values are standardised to equal None.
    """

    amount_factor_dict = {"ug": 0.000001, "mg": 0.001, "g": 1}
    mass_factor_dict = {"umol": 0.000001, "mmol": 0.001, "mol": 1}
    none_dict = {"-select-": None}

    def __init__(self, db_reaction: models.Reaction):
        """db_reaction is the reaction we are generating metadata for"""
        self.db_reaction = db_reaction
        self.rxn_data = json.loads(db_reaction.reaction_table_data)
        self.summary_data = json.loads(db_reaction.summary_table_data)
        # blank dictionary if 'to_export' is not present or convert the 'to_export' to standard dict if it is present
        self.exported_from_js = (
            self.list_to_dict(json.loads(self.summary_data["to_export"]))
            if self.summary_data.get("to_export")
            else {}
        )
        self.metadata_access_methods = self.get_metadata_access_methods()
        self.units = {}
        self.metadata_dict = {}
        self.make_metadata_dict()

    def make_metadata_dict(self):
        """Called on __init__, makes the metadata dictionary by iterating through the access method dict"""
        self.initialise_units()
        for key, access_method in self.metadata_access_methods.items():
            self.metadata_dict.update({key: access_method()})

    def get_dict(self) -> Dict:
        """Returns the reaction's metadata as a dictionary"""
        return self.metadata_dict

    # 'who' did the reaction

    def read_creator_username(self):
        """Returns the username of the reaction creator"""
        return self.db_reaction.creator_person.user.username

    def read_workbook(self):
        """Returns the name of the workbook the reaction was made in"""
        return self.db_reaction.workbook.name

    def read_workgroup(self):
        """Returns the name of the workgroup the reaction was made in"""
        return self.db_reaction.workbook.WorkGroup.name

    # 'when' was the reaction done

    def read_time_of_creation(self):
        """Returns the time the reaction was made"""
        return self.db_reaction.time_of_creation.isoformat()

    def read_reaction_completed(self):
        """Returns whether the reaction is 'complete' or 'not complete'"""
        return self.db_reaction.complete

    # what was the reaction
    def read_reaction_name(self) -> str:
        """Returns the user created reaction name"""
        return self.db_reaction.name

    def read_reaction_id(self) -> str:
        """Returns the autogenerated reaction id. e.g., RXN-001"""
        return self.db_reaction.reaction_id

    def read_solvents(self) -> Optional[List[Dict]]:
        """
        If solvents are present, returns a list of dictionaries with the solvent name, SMILES, volume and volume unit
        """
        if not self.check_compounds_present("solvent"):
            return None

        solvent_data = self.get_compound_data("solvent")
        self.standardise_compound_data(solvent_data)
        number_of_solvents = len(solvent_data["names"])
        solvent_data.update(
            {
                "volume_unit": number_of_solvents * [self.units["solvent_volume"]],
            }
        )
        return self.format_solvent_data(solvent_data)

    def read_reactants(self) -> Optional[List[Dict]]:
        """
        If reactants are present, returns a list of dictionaries with reactant name, SMILES, masses, amounts and units
        """
        if not self.check_compounds_present("reactants"):
            return None

        reactant_data = self.get_compound_data("reactant")
        self.standardise_compound_data(reactant_data)
        number_of_reactants = len(reactant_data["names"])
        reactant_data.update(
            {
                "mass_unit": number_of_reactants * [self.units["reactant_mass"]],
                "amount_unit": number_of_reactants * [self.units["reactant_amount"]],
            }
        )
        return self.format_compound_data(reactant_data)

    def read_reagents(self) -> Optional[List[Dict]]:
        """
        If reagents are present, returns a list of dictionaries with reagent name, SMILES, masses, amounts, and units
        """
        if not self.check_compounds_present("reagents"):
            return None
        reagent_data = self.get_compound_data("reagent")
        self.standardise_compound_data(reagent_data)
        number_of_reagents = len(reagent_data["names"])
        reagent_data.update(
            {
                "mass_unit": number_of_reagents * [self.units["reactant_mass"]],
                "amount_unit": number_of_reagents * [self.units["reactant_amount"]],
            }
        )
        return self.format_compound_data(reagent_data)

    def read_products(self) -> Optional[List[Dict]]:
        """
        If products are present, returns a list of dictionaries with product name, SMILES, masses, amounts, and units
        """
        if not self.check_compounds_present("products"):
            return None
        product_data = self.get_compound_data("product")
        self.standardise_compound_data(product_data)
        number_of_products = len(product_data["names"])
        product_data.update(
            {
                "mass_unit": number_of_products * [self.units["product_mass"]],
                "amount_unit": number_of_products * [self.units["product_amount"]],
            }
        )
        return self.format_compound_data(product_data)

    def read_experimental_writeup(self) -> str:
        """Returns the experimental writeup, currently overlaps with the description."""
        description = self.db_reaction.description
        return description if description != "" else None

    def read_standard_protocols_used(self) -> Optional[List[str]]:
        """Returns a list of the standard protocol fields that the user ticked that apply to this reaction"""
        standard_protocols = [
            ReactionStringMapper.radio_buttons(
                services.utils.camelCase_to_snake_case(x)
            )
            for x in self.summary_data.get("radio_buttons")
            if x is not None
        ]
        return standard_protocols if standard_protocols != [] else None

    def read_file_attachment_names(self) -> Optional[List[str]]:
        """Returns a list of the names for the file attachments for this reaction"""
        file_names = [x.display_name for x in self.db_reaction.file_attachments]
        return file_names if file_names != [] else None

    def read_addenda(self) -> Optional[List[Dict]]:
        """Returns the addenda a list of dictionaries with the author, comment text, and date."""
        return (
            [x.to_dict() for x in self.db_reaction.addenda]
            if self.db_reaction.addenda
            else None
        )

    # how sustainable was the reaction
    def read_solvent_sustainability(self) -> Optional[List[str]]:
        """Returns a list of strings to describe the sustainability of each solvent using the CHEM21 system"""
        cleaned_solvent_primary_keys = remove_default_data(self.db_reaction.solvent)
        sustainability_flags = [
            ReactionStringMapper.solvent_sustainability(
                services.solvent.sustainability_from_primary_key(x)
            )
            for x in cleaned_solvent_primary_keys
        ]
        return sustainability_flags if sustainability_flags != [] else None

    def read_reaction_safety(self) -> Optional[str]:
        """
        Returns a string to for overall risk of the reaction according to the hazard/exposure/risk matrix.
        L-low, M-Medium, H-High, VH-Very High
        """
        return self.exported_from_js.get("Overall-Risk-Rating")

    def read_temperature(self) -> Optional[str]:
        """Returns the reaction temperature the user entered into the summary table"""
        temperature = self.summary_data.get("reaction_temperature")
        return temperature if temperature != "" else None

    def read_element_sustainability(self) -> Optional[str]:
        """Returns the element sustainability expressed as a range of years"""
        return ReactionStringMapper.element_sustainability(
            self.summary_data.get("element_sustainability")
        )

    def read_batch_or_flow(self) -> Optional[str]:
        """Returns 'batch' or 'flow' from the user completed dropdown in the summary table"""
        batch_flow = self.summary_data.get("batch_flow")
        return batch_flow if batch_flow != "-select-" else None

    def read_purification_method(self) -> Optional[str]:
        """Returns the method of purification used from the dropdown the user selected in the summary table"""
        return ReactionStringMapper.isolation(self.summary_data.get("isolation_method"))

    def read_stoichiometry(self) -> Optional[str]:
        """Returns string to indicate whether catalyst was used or stoichiometric reagents or neither"""
        stoich = self.summary_data.get("catalyst_used")
        return stoich if stoich != "-select-" else None

    def read_catalyst_recovery(self) -> Optional[str]:
        """Returns string to indicate, if a catalyst was used, was it recovered."""
        cat_recovery = self.summary_data.get("catalyst_recovered")
        return cat_recovery if cat_recovery != "-select-" else None

    def read_atom_economy(self) -> Optional[str]:
        """
        Returns the atom economy as a string - percentage of mass of starting atoms in desired product
        Based on the theoretical yield hence does not require actual yield data
        """
        return (
            self.exported_from_js.get("Atom Efficiency")
            if self.exported_from_js.get("Atom Efficiency")
            else None
        )

    def read_mass_efficiency(self) -> str:
        """
        Returns the mass efficiency as a string - percentage of mass of starting atoms in desired product.
        Based on the real yield - does require user entered mass yield
        """
        return (
            self.summary_data.get("mass_efficiency")
            if self.summary_data.get("mass_efficiency")
            else None
        )

    def read_percent_yield(self) -> Optional[str]:
        """
        Returns the percentage yield for the reaction if the user has entered the product mass they obtained

        Returns:
              percentage yield as a string
        """
        actual_mass_yield = self.summary_data.get("real_product_mass")
        # a yield of zero is a valid yield
        if not actual_mass_yield and actual_mass_yield != "0":
            return None
        return str(
            services.sustainability.percent_yield_from_rxn_data(
                self.rxn_data, float(actual_mass_yield)
            )
        )

    def read_conversion(self) -> Optional[str]:
        """Returns the conversion from the metrics if present else None"""
        return (
            self.summary_data.get("conversion")
            if self.summary_data.get("conversion")
            else None
        )

    def read_selectivity(self):
        """Returns the selectivity from the metrics if present else None"""
        return (
            self.summary_data.get("selectivity")
            if self.summary_data.get("conversion")
            else None
        )

    def read_sustainability_flags(self) -> Dict:
        """Returns a dictionary with a colour-coded flag for each metric. Green, amber, red or None"""
        return {
            "temperature_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Temperature Sustainability")
            ),
            "element_sustainability_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Elements Sustainability")
            ),
            "batch_or_flow_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Batch or Flow Sustainability")
            ),
            "isolation_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Isolation Sustainability")
            ),
            "catalyst_used_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Catalyst Sustainability")
            ),
            "catalyst_recovery_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Recovery Sustainability")
            ),
            "atom_economy_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Atom Economy Sustainability")
            ),
            "mass_efficiency_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Mass Efficiency Sustainability")
            ),
            "yield_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Yield Sustainability")
            ),
            "conversion_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Conversion Sustainability")
            ),
            "selectivity_flag": ReactionStringMapper.css_classes_to_chem21_colour_flag(
                self.exported_from_js.get("Selectivity Sustainability")
            ),
        }

    def initialise_units(self):
        """Updates the units dictionary with all units used in the reaction."""
        self.units["reactant_mass"] = self.rxn_data["mass_units"]
        self.units["reactant_amount"] = self.rxn_data["amount_units"]
        self.units["product_amount"] = self.rxn_data["product_amount_units"]
        self.units["product_mass"] = self.rxn_data["product_mass_units"]
        self.units["solvent_volume"] = self.rxn_data["solvent_volume_units"]

    def check_compounds_present(self, component: str) -> bool:
        """
        component: reactant/reagent/solvent/product
        Returns True if that compound type is present in the reaction. e.g., A reaction has no solvents returns False
        """
        return bool(remove_default_data(getattr(self.db_reaction, component)))

    @staticmethod
    def standardise_compound_data(compound_data: Dict):
        """
        Sometimes it is possible for an empty string or other erroneous variables to enter the json
        In this function we trim the length of lists in the dictionary so all are equal to the shortest one
        """
        length_of_values = []
        for value in compound_data:
            length_of_values.append(len(value))

        length_of_values = [len(x) for x in compound_data.values()]
        if min(length_of_values) == max(length_of_values):
            return
        # if there are lists which are longer than others then we trim the list to be of the same length
        min_length = min(length_of_values)
        for key, value in compound_data.items():
            if len(value) > min_length:
                compound_data[key] = value[:min_length]

    @staticmethod
    def format_compound_data(compound_data: Dict) -> List[Dict]:
        """
        Formats compound data from a dictionary containing lists into a list of dictionaries. Not for solvents

        Args:
            compound_data - dict with all the values of equal length for compound data
        Returns:
             A list of dictionaries, 1 per compound for a specific reaction component.
        """

        return [
            {
                "name": name,
                "SMILES": smiles,
                "amount": amount,
                "mass": mass,
                "amount_unit": amount_unit,
                "mass_unit": mass_unit,
                "h_codes": h_code,
                "physical_form": physical_form,
            }
            for name, smiles, amount, mass, amount_unit, mass_unit, h_code, physical_form in zip(
                compound_data["names"],
                compound_data["smiles"],
                compound_data["amounts"],
                compound_data["masses"],
                compound_data["amount_unit"],
                compound_data["mass_unit"],
                compound_data["h_codes"],
                compound_data["physical_forms"],
            )
        ]

    @staticmethod
    def format_solvent_data(solvent_data: Dict) -> List[Dict]:
        """
        Formats solvent data from a dictionary containing lists into a list of dictionaries. Only for solvents

        Args:
            solvent_data - dict with all the values of equal length for solvent data
        Returns:
             A list of dictionaries, 1 per solvent.
        """
        return [
            {
                "name": name,
                "SMILES": smiles,
                "volume": volume,
                "volume_unit": volume_unit,
                "physical_form": physical_form,
                "hazard_codes": h_code,
            }
            for name, smiles, volume, volume_unit, physical_form, h_code in zip(
                solvent_data["names"],
                solvent_data["smiles"],
                solvent_data["volumes"],
                solvent_data["volume_unit"],
                solvent_data["physical_forms"],
                solvent_data["h_codes"],
            )
        ]

    def get_compound_data(self, component: str) -> Dict:
        """
        For a component reagent/solvent/reactant/product we make a dictionary with the required data
        Args:
            component - indicate what reaction role we are making the dict for (e.g., reactant)
        Returns:
            A dictionary with compound data

        """
        compound_dict = {
            "names": self.get_compound_names(component),
            "h_codes": self.get_h_codes(component),
            "physical_forms": self.get_physical_forms(component),
        }
        if component == "solvent":
            compound_dict["volumes"] = self.get_volumes(component)
            compound_dict["smiles"] = services.all_compounds.get_smiles_list(
                self.db_reaction.solvent
            )
        else:  # for reactants/reagents/products
            compound_dict["smiles"] = self.get_compound_smiles(
                component + "s"
            )  # plural used in models.reaction.components
            compound_dict["masses"] = self.get_masses(component)
            compound_dict["amounts"] = self.get_amounts(component)
        return compound_dict

    def get_compound_names(self, component: str) -> List[str]:
        """Returns a list of compound names for type of component: reactant/reagent/solvent/product"""
        return [x for x in self.rxn_data[f"{component}_names"]]

    def get_compound_smiles(self, component: str) -> List[str]:
        """Returns a list of SMILES for type of component: reactant/reagent/product. Not for Solvents"""
        return getattr(self.db_reaction, component)

    def get_amounts(self, component: str) -> List[str]:
        """Returns a list of amounts for type of component: reactant/reagent/solvent/product"""
        return [
            str(round(float(x), 2)) for x in self.rxn_data[f"{component}_amounts_raw"]
        ]

    def get_masses(self, component: str) -> List[str]:
        """Returns a list of masses for type of component: reactant/reagent/solvent/product"""
        return [
            str(round(float(x), 2)) for x in self.rxn_data[f"{component}_masses_raw"]
        ]

    def get_volumes(self, component) -> List[str]:
        """Returns a list of volumes for type of component: reactant/reagent/solvent/product"""
        return [str(round(float(x), 2)) for x in self.rxn_data[f"{component}_volumes"]]

    def get_h_codes(self, component) -> List[str]:
        """Returns a list of GHS hazard codes for type of component: reactant/reagent/solvent/product"""
        return [x for x in self.rxn_data[f"{component}_hazards"]]

    def get_physical_forms(self, component):
        """Returns a list of physical forms for type of component: reactant/reagent/solvent/product"""
        return [
            ReactionStringMapper.physical_forms(x)
            for x in self.rxn_data[f"{component}_physical_forms"]
        ]

    def get_metadata_access_methods(self) -> Dict:
        """Returns a dictionary with the key and the associated access function used to get the value"""
        return {
            # who
            "creator_username": self.read_creator_username,
            "creator_workbook": self.read_workbook,
            "creator_workgroup": self.read_workgroup,
            # when
            "time_of_creation": self.read_time_of_creation,
            "reaction_completed": self.read_reaction_completed,
            # what
            "reaction_id": self.read_reaction_id,
            "reaction_name": self.read_reaction_name,
            "reagents": self.read_reagents,
            "reactants": self.read_reactants,
            "solvents": self.read_solvents,
            "products": self.read_products,
            "experimental_writeup": self.read_experimental_writeup,
            "standard_protocols_used": self.read_standard_protocols_used,
            "file_attachment_names": self.read_file_attachment_names,
            "addenda": self.read_addenda,
            # sustainability
            "solvent_sustainability": self.read_solvent_sustainability,
            "safety": self.read_reaction_safety,
            "temperature": self.read_temperature,
            "element_sustainability": self.read_element_sustainability,
            "batch_or_flow": self.read_batch_or_flow,
            "purification_method": self.read_purification_method,
            "catalyst_use": self.read_stoichiometry,
            "catalyst_recovery": self.read_catalyst_recovery,
            "atom_economy": self.read_atom_economy,
            "mass_efficiency": self.read_mass_efficiency,
            "percent_yield": self.read_percent_yield,
            "conversion": self.read_conversion,
            "selectivity": self.read_selectivity,
            "sustainability_data": self.read_sustainability_flags
            # [{'batch_or_flow': {'value': 'batch', 'flag': 1 or 'green'}]
            # Add other metadata properties with corresponding access methods...
        }

    @staticmethod
    def list_to_dict(ls: List[Dict]) -> Dict:
        """From [{'key': 'risk-score', 'value': 'danger'}] -> {'risk-score': 'danger'}"""
        result_dict = {}
        for item in ls:
            result_dict[item["key"]] = item["value"]
        return result_dict


class ReactionStringMapper:
    @staticmethod
    def chem21_to_numerical(chem21_str: str) -> str:
        """Used to replace the css classes that indicate colour with a numerical equivalent"""
        chem21_dict = {
            "recommended": "1",
            "problematic": "2",
            "hazard-warning": "3",
            "hazard-hazardous": "4",
        }
        return chem21_dict[chem21_str]

    @staticmethod
    def css_classes_to_chem21_colour_flag(css_class: str) -> str:
        css_chem21_dict = {
            "hazard-reset-hazard": None,
            "hazard-acceptable": "green",
            "hazard-warning": "amber",
            "hazard-hazardous": "red",
            None: None,
        }
        return css_chem21_dict[css_class]

    @staticmethod
    def isolation(isolation_method: str) -> str:
        """Uses to replace the isolation dropdown index with the corresponding string value"""
        isolation_dict = {
            "0": "",
            "1": "Column",
            "2": "HPLC",
            "3": "Ion exchange",
            "4": "Crystallization",
            "5": "Filtration",
            "6": "Multiple recryst.",
            "7": "Distillation < 140 degC",
            "8": "Distillation > 140 degC",
            "undefined": None,
        }
        return isolation_dict[isolation_method]

    @staticmethod
    def element_sustainability(element_value: str) -> str:
        """Replaces the element sustainability dropdown index with the corresponding string value"""
        elements_dict = {
            "3": "+500 years",
            "2": "50-500 years",
            "1": "5-50 years",
            "0": "",
            "undefined": None,
        }
        return elements_dict[element_value]

    @staticmethod
    def solvent_sustainability(solvent_flag: int) -> str:
        """Replaces the solvent sustainability numerical flag found in the database with the corresponding string"""
        # get solvent sustainability
        solvent_flag_dict = {
            4: "Recommended",
            3: "Problematic",
            2: "Hazardous",
            1: "Highly Hazardous",
            5: "Unknown",  # non-chem21 elsewhere in code sometimes.
            "undefined": None,
        }
        return solvent_flag_dict[solvent_flag]

    @staticmethod
    def physical_forms(physical_form):
        """Replaces the physical form dropdown index with the corresponding string value"""
        physical_form_dict = {
            "0": "",
            "1": "Dense Solid",
            "2": "Non-volatile liquid (b.p. > 130 degC)",
            "3": "Unknown",
            "4": "Dusty Solid",
            "5": "Lyophilised solid",
            "6": "Volatile liquid (b.p. 70-130 degC)",
            "7": "Gas",
            "8": "Highly volatile liquid (b.p. < 70 degC)",
            "9": "Aerosol",
            "10": "Solution that promotes skin absorption",
            "undefined": None,
        }
        return physical_form_dict[physical_form]

    @staticmethod
    def radio_buttons(radio_button) -> str:
        """
        We take the radio button name and change it if it appears in the dict, otherwise we return the original

        Args:
            radio_button - the string from the frontend of the selected radio button

        Returns:
            Either the new string from the dict value or the original if the string was not in the dictionary,
        """
        radio_button_dict = {
            "slight": "hazard: slight",
            "serious": "hazard: serious",
            "major": "hazard: major",
            "low_likelihood": "risk: low likelihood",
            "possible": "risk: possible",
            "frequent_occur": "risk: frequent occur",
            "individual": "consequences: individual",
            "local_labs": "consequences: local labs",
            "building_wide": "consequences: building wide",
        }
        new_string = radio_button_dict.get(radio_button)
        if new_string:
            return new_string
        return radio_button


def get_metadata_keys() -> List[str]:
    """
    Returns the keys used in the metadata
    """
    return [
        "temperature",
        "solvents",
        "reagents",
        "creator_username",
        "creator_workbook",
        "creator_workgroup",
        "time_of_creation",
        "reaction_completed",
        "purification_method",
        "batch_or_flow",
        "percent_yield",
        "reactants",
        "experimental_writeup",
        "standard_protocols_used",
        "nmr_data",
        "file_attachment_names",
    ]


def remove_default_data(list_: List):
    """Removes default database values from lists for database reaction components"""
    return [x for x in list_ if x not in ["{", "}", ""]]
