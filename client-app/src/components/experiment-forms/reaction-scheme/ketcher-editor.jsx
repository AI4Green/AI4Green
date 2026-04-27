import { Button, FormControl, HStack, Text, VStack } from "@chakra-ui/react";
import { FormHelpError } from "components/core/forms";
import { useBackendApi } from "contexts";
import { useField, useFormikContext } from "formik";
import {
  reactionSmilesToReactantsAndProductsSmiles,
  replaceSmilesSymbols,
} from "helpers";
import { useRef, useState } from "react";
import { FaSync } from "react-icons/fa";
import { GiMaterialsScience } from "react-icons/gi";

const KETCHER_IFRAME_SRC = "/static/ketcher/index.html";
const KETCHER_EDITOR_INITALS_VALUES = {
  sketcherSmiles: "", // smiles from the Ketcher
  reactants: [], // reactants extracted using the smiles
  products: [], // products extracted using the smiles
  smiles: "", // smiles from the Ketcher with replaced symbols
  data: null, // response from the reaction table api
};

export const KetcherEditor = ({ parentName, name, isDisabled }) => {
  const [field, , helpers] = useField(name);
  const { setFieldValue } = useFormikContext();

  const [isLoading, setIsLoading] = useState();
  const [feedback, setFeedback] = useState();

  const ketcherIframe = useRef(null);

  const { reactionTable: action } = useBackendApi();

  const ketcherWindow = ketcherIframe.current?.contentWindow;

  const handleKetcherOnLoad = async () => {
    const { sketcherSmiles } = field.value || KETCHER_EDITOR_INITALS_VALUES;

    // keep checking until the ketcher is loaded
    while (!ketcherWindow.ketcher) {
      await new Promise((resolve) => setTimeout(resolve, 100)); // wait 100ms before checking again
    }

    if (sketcherSmiles) {
      try {
        await ketcherWindow.ketcher.setMolecule(sketcherSmiles);
      } catch (error) {
        console.error("Error setting molecule in Ketcher editor:", error);
      }
    }
  };

  // Extract smiles from sketcher and use it to get data from the reaction table data
  const handleDataGenerate = async () => {
    if (!ketcherWindow) return;

    const sketcherSmiles = await ketcherWindow.ketcher.getSmiles();
    const { reactants, products } =
      reactionSmilesToReactantsAndProductsSmiles(sketcherSmiles);

    const smiles = replaceSmilesSymbols(sketcherSmiles);

    if (!reactants || !products || !smiles) {
      setFeedback("Invalid reaction provided.");
      helpers.setValue(KETCHER_EDITOR_INITALS_VALUES);
      return;
    }

    try {
      setIsLoading(true);
      const data = await action.process(reactants, products, smiles);
      const reactionImage = await ketcherWindow.ketcher.generateImage(
        sketcherSmiles,
      );

      helpers.setValue({
        sketcherSmiles,
        reactionImage: {
          ...field.value?.reactionImage,
          image: reactionImage,
          isNew: !field.value?.reactionImage,
          isMarkedForDeletion: false,
        },
        reactants,
        products,
        smiles,
        data,
      });
      feedback && setFeedback(null);
    } catch (error) {
      setFeedback(error?.message ?? "Something went wrong.");
      helpers.setValue(KETCHER_EDITOR_INITALS_VALUES);
    } finally {
      setIsLoading(false);
    }
  };

  // reset the ketcher and the reaction table values also the ketcher editor values
  const handleKetcherReset = async () => {
    if (ketcherWindow) {
      await ketcherWindow.ketcher.setMolecule("");
    }
    helpers.setValue({
      ...KETCHER_EDITOR_INITALS_VALUES,
      reactionImage: field.value?.reactionImage
        ? {
            ...field.value?.reactionImage,
            isMarkedForDeletion: true,
          }
        : null,
    });
    setFieldValue(`${parentName}.reactionTable`, []);
  };

  return (
    <FormControl isRequired id={field.name} isInvalid={Boolean(feedback)}>
      <VStack minW="full" align="flex-start">
        <Text as="b">Reaction Sketcher</Text>
        <iframe
          ref={ketcherIframe}
          src={KETCHER_IFRAME_SRC}
          title="Ketcher App"
          width="100%"
          height="500px"
          allowFullScreen
          onLoad={handleKetcherOnLoad}
          style={{ pointerEvents: isLoading || isDisabled ? "none" : "auto" }}
        />
        <HStack>
          {!isDisabled && !field.value?.sketcherSmiles && (
            <Button
              leftIcon={<GiMaterialsScience />}
              isLoading={isLoading}
              disabled={isLoading}
              colorScheme="purple"
              size="sm"
              onClick={handleDataGenerate}
            >
              Generate reaction data
            </Button>
          )}
          {!isDisabled && field.value?.sketcherSmiles && (
            <Button
              leftIcon={<FaSync />}
              isLoading={isLoading}
              disabled={isLoading}
              colorScheme="orange"
              size="sm"
              onClick={handleKetcherReset}
            >
              Reset
            </Button>
          )}
        </HStack>
        <FormHelpError
          isInvalid={Boolean(feedback)}
          error={feedback}
          collapseEmpty
          replaceHelpWithError
        />
      </VStack>
    </FormControl>
  );
};
