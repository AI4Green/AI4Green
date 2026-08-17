import {
  Box,
  Button,
  FormControl,
  HStack,
  Image,
  Text,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { FormHelpError, FormikInput } from "components/core/forms";
import { LoadingIndicator } from "components/core/loading-indicator";
import { Modal } from "components/core/modal";
import { useBackendApi } from "contexts";
import { Form, Formik, useField } from "formik";
import { ContentPage } from "pages/content";
import { useState } from "react";
import { LuLoaderPinwheel } from "react-icons/lu";

export const ReactionPredictor = ({
  name,
  label = "Input SMILES",
  isDisabled,
}) => {
  const [isLoading, setIsLoading] = useState();
  const [feedback, setFeedback] = useState();
  const { isOpen, onOpen, onClose } = useDisclosure();

  const [field, , helpers] = useField(name);
  const { predictions: action } = useBackendApi();

  const handlePrediction = async ({ smiles }) => {
    setIsLoading(true);
    try {
      const result = await action.forwardPrediction(smiles);
      const data = await result.json();

      helpers.setValue({ smiles, predictions: data });
      feedback && setFeedback(null);
    } catch (e) {
      const error = await e.response?.json();
      setFeedback(error?.message || "Something went wrong.");
    } finally {
      setIsLoading(false);
    }
  };

  const modalBody = <ContentPage contentKey="reactionprediction" />;

  return (
    <Box w="full" align="flex-start">
      <FormControl id={field.name} isInvalid={Boolean(feedback)}>
        <Formik
          enableReinitialize
          initialValues={{
            ...(field.value || { smiles: "", predictions: [] }),
          }}
          onSubmit={handlePrediction}
        >
          {({ handleSubmit, values }) => (
            <Box w="full" borderRadius={7} borderWidth={1} p={4}>
              <Form>
                <VStack align="flex-start" spacing={8}>
                  <Button
                    size="xs"
                    onClick={onOpen}
                    leftIcon={<LuLoaderPinwheel />}
                    variant="link"
                  >
                    Learn about Reaction Predictions
                  </Button>
                  {isOpen && (
                    <Modal
                      body={modalBody}
                      title="Learn about Reaction Predictions"
                      onAction={onClose}
                      isOpen={isOpen}
                      contentMaxW="80vw"
                      contentMaxH="90vh"
                      bodyMaxH="70vh"
                      bodyOverflowY="auto"
                      cancelBtnEnable={false}
                    />
                  )}
                  <FormikInput
                    name="smiles"
                    label={label}
                    isDisabled={isDisabled}
                    isRequired
                  />
                  <Button
                    isLoading={isLoading}
                    loadingText="Predicting"
                    onClick={handleSubmit}
                    colorScheme="green"
                    size="sm"
                    px={4}
                    leftIcon={<LuLoaderPinwheel />}
                    isDisabled={values.smiles === "" || isLoading || isDisabled}
                    hidden={isDisabled}
                  >
                    Predict
                  </Button>
                  <FormHelpError
                    isInvalid={Boolean(feedback)}
                    error={feedback}
                    collapseEmpty
                    replaceHelpWithError
                  />
                </VStack>
              </Form>

              {field.value.predictions && (
                <Product
                  products={field.value.predictions}
                  isLoading={isLoading}
                />
              )}
            </Box>
          )}
        </Formik>
      </FormControl>
    </Box>
  );
};

const Product = ({ products, isLoading }) => {
  if (isLoading) {
    return <LoadingIndicator verb="Predicting" noun="products" />;
  }

  return (
    <VStack align="flex-start" py={8}>
      <Text fontWeight="semibold" fontSize="xl">
        Top 5 Predicted Products
      </Text>
      {products.map((product) => (
        <VStack
          key={product.id}
          align="flex-start"
          borderBottomWidth={1}
          borderRadius={5}
          p={4}
          fontSize="sm"
        >
          <HStack>
            <Text fontWeight="semibold">IUPAC Name:</Text>
            <Text>{product.iupacName}</Text>
          </HStack>
          <HStack>
            <Text fontWeight="semibold">Product SMILES:</Text>
            <Text>{product.product}</Text>
          </HStack>
          <HStack>
            <Text fontWeight="semibold">Score:</Text>
            <Text>{product.score}</Text>
          </HStack>
          <Image
            src={`data:image/png;base64,${product.reactionImage}`}
            alt={product.iupacName}
          />
          <HStack>
            <Text fontWeight="semibold">Synonyms:</Text>
            <Text>{product.synonyms.slice(0, 3).join(", ")}</Text>
          </HStack>
        </VStack>
      ))}
    </VStack>
  );
};
