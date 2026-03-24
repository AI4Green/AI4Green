export const replaceSmilesSymbols = (reaction) => {
  reaction = reaction.replace(/\+/g, "plus");
  reaction = reaction.replace(/-/g, "minus");
  return reaction.replace(/#/g, "sharp");
};

export const reactionSmilesToReactantsAndProductsSmiles = (sketcherSmiles) => {
  let reaction = sketcherSmiles.split(" |")[0];
  reaction = replaceSmilesSymbols(reaction);
  let array = reaction.split(">");
  let reactants = array[0]?.split(".");
  let products = array[2]?.split(".");
  return { reactants, products };
};

export const removeReagentsFromSmiles = (sketcherSmiles) => {
  let smiles2 = sketcherSmiles.split(">").slice(-1);
  let smiles1 = sketcherSmiles.split(">")[0];
  return smiles1 + ">>" + smiles2;
};
