import { Box, HStack, Text, useToast, VStack } from "@chakra-ui/react";
import { GLOBAL_PARAMETERS } from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useState } from "react";
import { object, string } from "yup";
import { TextField } from "./form-fields.jsx";
import { ErrorMessage, INV, InvButton, InvModal, Muted } from "./ui.jsx";
import { describeError } from "./use-inventory-action.js";

const FORM_ID = "inv-cheminventory-form";
const MAX_TOKEN_LENGTH = 512;

const schema = object({
  token: string()
    .trim()
    .required("Enter an API token")
    .max(MAX_TOKEN_LENGTH, "That doesn't look like an API token")
    .matches(/^\S+$/, "API tokens don't contain spaces"),
});

/**
 * Connect a workgroup to ChemInventory with an API token, replace the token,
 * or disconnect. The token is sent once and never shown again. The backend
 * only returns its last few characters.
 * @param {object} props
 * @param {string} props.workgroupName The workgroup.
 * @param {object} props.connection The connection status.
 * @param {(isConnected: boolean) => void} props.onChanged Called after
 *   connecting or disconnecting.
 * @param {() => void} props.onClose Closes the modal.
 */
export const ChemInventoryModal = ({
  workgroupName,
  connection,
  onChanged,
  onClose,
}) => {
  const { inventory } = useBackendApi();
  const toast = useToast();
  const [isLoading, setIsLoading] = useState(false);
  const [error, setError] = useState(null);
  const [isConfirmingDisconnect, setIsConfirmingDisconnect] = useState(false);
  const isConnected = Boolean(connection.configured);

  /**
   * Runs a request, showing a toast if it works and an error if it doesn't.
   * @returns {Promise<boolean>} Whether it worked.
   */
  const run = async (request, successTitle) => {
    setIsLoading(true);
    setError(null);
    try {
      await request();
      toast({
        title: successTitle,
        status: "success",
        duration: GLOBAL_PARAMETERS.ToastDuration,
        isClosable: true,
        position: "top",
      });
      return true;
    } catch (e) {
      // not logged: the error holds the request options, including the token
      setError(await describeError(e));
      return false;
    } finally {
      setIsLoading(false);
    }
  };

  const handleSubmit = async ({ token }) => {
    const ok = await run(
      () => inventory.saveChemInventoryConnection(workgroupName, token.trim()),
      "ChemInventory connected",
    );
    if (ok) {
      onChanged(true);
      onClose();
    }
  };

  const handleDisconnect = async () => {
    const ok = await run(
      () => inventory.deleteChemInventoryConnection(workgroupName),
      "ChemInventory disconnected",
    );
    if (ok) {
      onChanged(false);
      onClose();
    }
  };

  return (
    <InvModal
      isOpen
      onClose={onClose}
      title="ChemInventory"
      actionLabel={isConnected ? "Replace token" : "Connect"}
      formId={FORM_ID}
      isLoading={isLoading}
    >
      {error && <ErrorMessage>{error}</ErrorMessage>}

      {isConnected && (
        <Box
          mb="16px"
          p="12px"
          borderRadius="8px"
          border="1px solid"
          borderColor={INV.line}
          bg={INV.bg}
        >
          <Text m={0}>
            Connected
            {connection.cheminventoryUser &&
              ` as ${connection.cheminventoryUser}`}
            {connection.inventoryName && ` to ${connection.inventoryName}`}.
          </Text>
          <Muted mt="4px">
            Token ending {connection.tokenHint}
            {connection.connectedBy && `, added by ${connection.connectedBy}`}
          </Muted>
          {isConfirmingDisconnect ? (
            <Box mt="12px">
              <Text m={0} mb="8px">
                Disconnect ChemInventory? The stored token will be deleted, and
                nobody in the workgroup will see ChemInventory stock here until
                a new one is added.
              </Text>
              <HStack spacing="8px">
                <InvButton
                  onClick={() => setIsConfirmingDisconnect(false)}
                  isDisabled={isLoading}
                >
                  Keep connected
                </InvButton>
                <InvButton
                  danger
                  onClick={handleDisconnect}
                  isLoading={isLoading}
                >
                  Disconnect
                </InvButton>
              </HStack>
            </Box>
          ) : (
            <InvButton
              mt="12px"
              onClick={() => setIsConfirmingDisconnect(true)}
              isDisabled={isLoading}
            >
              Disconnect…
            </InvButton>
          )}
        </Box>
      )}

      <Formik
        initialValues={{ token: "" }}
        validationSchema={schema}
        onSubmit={handleSubmit}
      >
        <VStack
          as={Form}
          id={FORM_ID}
          noValidate
          spacing="16px"
          align="stretch"
        >
          <TextField
            name="token"
            label={isConnected ? "New API token" : "API token"}
            type="password"
            // browsers ignore "off" on password fields; "new-password" stops a
            // password manager filling in the user's AI4Green password, which
            // would then be sent to ChemInventory
            autoComplete="new-password"
            spellCheck={false}
            isRequired
            autoFocus={!isConnected}
            helpText="Everyone in this workgroup will see the stock this ChemInventory user can see."
          />
          <Muted fontSize="13px">
            A ChemInventory Group Administrator can create a token under
            Inventory Management → Inventory Tools → API Tokens. Create it for a
            read-only ChemInventory user rather than an administrator, so it
            can't change anything in ChemInventory. The token is stored
            encrypted and can't be viewed again once saved.
          </Muted>
        </VStack>
      </Formik>
    </InvModal>
  );
};
