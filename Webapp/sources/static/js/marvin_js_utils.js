function newMarvinSketcher(IDSelector) {
    return ChemicalizeMarvinJs.createEditor(IDSelector)
}

function exportImageString(marvin) {
    // saves image to a hidden html input
    let settings = {'width': 600, 'height': 400};
    marvin.exportStructure("jpeg", settings).then(function (source) {
        $("#js-reloaded-image").val(source);
    });
}


function replaceSmilesSymbols(reaction){
    //replace the plus, minus and sharp signs from a smiles string to trasfer it properly
    reaction = reaction.replace(/\+/g, 'plus');
    reaction = reaction.replace(/-/g, 'minus');
    return reaction.replace(/#/g, 'sharp');
}

function processSmiles(smiles) {
    let reaction = smiles.split(" |")[0];
    reaction = replaceSmilesSymbols(reaction)
    //split into reactants, solvents and product
    let array = reaction.split(">");
    let reactants = array[0].split(".");
    let rl = reactants.length;
    let products = array[2].split(".");
    let reagents = array[1].split(".");
    let compounds, rgl, ending, group, element, c0, c, r, rg;
    let i, k, p;
    if (reagents[0] !== "") {
        rgl = reagents.length;
        compounds = reactants.concat(reagents, products);
    } else {
        rgl = 0;
        compounds = reactants.concat(products);
    }
    // if "f:" is at the end of the smiles string, it contains an ionic compound
    if (smiles.includes("f:") === true) {
        [reactants, products] = translateSalts(smiles, reactants, products)
        } else {
            ending = "none";
    }
    return [reactants, products]
}

function translateSalts(smiles, reactants, products){
    // split after the comma to see how many ionic compounds in reaction smiles
    let ionList = smiles.split("f:")[1].split('|')[0].split(',');
    // for each ionic compound, join all ions to make the full salt and replace ions in reactant/product list
    let productSaltDictList = []
    let reactantSaltDictList = []
    for (let adjacentIons of ionList) {
        // make list of index, and smiles, then join to make a string of the full salt from the ions
        let ionIndexList = adjacentIons.split('.').map(Number);
        let ionCompoundsSmiles = [] // [Pd+, OAc-, OAc-]
        for (let idx of ionIndexList) {
            ionCompoundsSmiles.push(compounds[idx])
        }
        let saltSmiles = ionCompoundsSmiles.join(".")
        // determine if reactant or product
        let rxnComponent
        if (rl > ionIndexList[0]) {
            rxnComponent = "reactant"
        } else {
            rxnComponent = "product"
        }
        // ions in previous salts for only same rxn component (reactant or product) - to determine insert position
        // when splicing into reactant/product list.
        let ionsInPreviousSalts = 0
        let insertPosition
        if (rxnComponent === 'reactant') {
            for (let previousSaltDict of reactantSaltDictList) {
                ionsInPreviousSalts += previousSaltDict["numberOfIons"]
            }
            insertPosition = ionIndexList[0] - ionsInPreviousSalts + reactantSaltDictList.length
            reactantSaltDictList.push({
                "saltSmiles": saltSmiles,
                "numberOfIons": ionCompoundsSmiles.length,
                "insertPosition": insertPosition
            })
        }
        if (rxnComponent === "product") {
            for (let previousSaltDict of productSaltDictList) {
                ionsInPreviousSalts += previousSaltDict["numberOfIons"]
            }
            insertPosition = ionIndexList[0] - ionsInPreviousSalts + productSaltDictList.length - rl
            productSaltDictList.push({
                "saltSmiles": saltSmiles,
                "numberOfIons": ionCompoundsSmiles.length,
                "insertPosition": insertPosition
            })
        }
    }
    // splice the salts into the lists replacing the ions that they contain
    for (let saltDict of reactantSaltDictList) {
        reactants.splice(saltDict["insertPosition"], saltDict["numberOfIons"], saltDict["saltSmiles"])
    }
    for (let saltDict of productSaltDictList){
            products.splice(saltDict["insertPosition"], saltDict["numberOfIons"], saltDict["saltSmiles"])
    }
    return [reactants, products]
}

function removeReagentsFromSmiles(smiles){
    // remove reagents from reaction smiles
    let smiles2 = smiles.split('>').slice(-1);
    let smiles1 = smiles.split('>')[0];
    return smiles1 + '>>' + smiles2;
}

function exportSmilesFromSketcher(marvin) {
    return marvin.exportStructure('cxsmiles', {'extra': 'f'}).then(function (smiles) {
        return smiles
    });
}