import {
  Badge as ChakraBadge,
  Box,
  Button,
  Flex,
  Heading,
  HStack,
  IconButton,
  Modal,
  ModalBody,
  ModalCloseButton,
  ModalContent,
  ModalFooter,
  ModalHeader,
  ModalOverlay,
  Td,
  Text,
  Th,
} from "@chakra-ui/react";

/*
 * Shared styling for the inventory page. Colours and sizes are taken from the
 * original inventory.css mockup. Variants and sizes are always set explicitly
 * so the page doesn't shift if the app theme's defaults change.
 */

export const INV = {
  accent: "#17543d",
  accentHover: "#0f3d2c",
  accentSoft: "#e3f1ea",
  ink: "#1b2a24",
  muted: "#56665e",
  line: "#dce3de",
  lineStrong: "#c4cfc8",
  // borders of interactive controls; the lightest shade in this palette with
  // 3:1 contrast against every background they sit on (WCAG 1.4.11)
  control: "#7e8c84",
  rowLine: "#e8ede9",
  bg: "#f6f7f4",
  rowOpen: "#fafbf9",
  sub: "#f1f5f2",
  subLine: "#e1e8e3",
  danger: "#8a2a1b",
  dangerHover: "#6e2115",
};

// replaces Chakra's box-shadow focus ring with the mockup's outline
export const focusRing = {
  outline: "2px solid",
  outlineColor: INV.accent,
  outlineOffset: "2px",
  boxShadow: "none",
};

// marks focused, hovered and selected states without relying on colour (WCAG 1.4.1)
export const underline = {
  textDecoration: "underline",
  textDecorationThickness: "2px",
  textUnderlineOffset: "4px",
};

// matches the mockup's (max-width: 720px) breakpoint exactly
export const MOBILE = "@media (max-width: 720px)";

// variant="unstyled" still gets Chakra's base and size styles, so reset those
const bareButton = {
  variant: "unstyled",
  display: "inline-flex",
  h: "auto",
  minW: 0,
  minH: "44px",
  lineHeight: "inherit",
  _focusVisible: { ...focusRing, ...underline },
};

// shared look for text inputs, selects and textareas
export const controlProps = {
  variant: "outline",
  h: "44px",
  fontSize: "14px",
  borderRadius: "8px",
  borderColor: INV.lineStrong,
  bg: "white",
  color: INV.ink,
  focusBorderColor: INV.accent,
  _hover: { borderColor: INV.lineStrong },
};

// solid colours for the filled button styles, keyed by InvButton's tone
const SOLID = {
  primary: { bg: INV.accent, hover: INV.accentHover },
  danger: { bg: INV.danger, hover: INV.dangerHover },
};

/**
 * Outlined button, or a solid button with `primary` (accent) or `danger` (red).
 * @param {object} props
 * @param {boolean} [props.primary] Use the solid accent style.
 * @param {boolean} [props.danger] Use the solid red style, for destructive actions.
 */
export const InvButton = ({ primary = false, danger = false, ...props }) => {
  const solid = danger ? SOLID.danger : primary ? SOLID.primary : null;
  const bg = solid ? solid.bg : "white";
  const borderColor = solid ? solid.bg : INV.control;
  return (
    <Button
      {...bareButton}
      px="16px"
      borderRadius="8px"
      border="1px solid"
      borderColor={borderColor}
      bg={bg}
      color={solid ? "white" : INV.ink}
      fontSize="14px"
      fontWeight={solid ? 500 : 400}
      _disabled={{ opacity: 0.5, cursor: "default" }}
      _hover={{
        bg: solid ? solid.hover : INV.sub,
        borderColor: solid ? solid.hover : INV.accent,
        ...underline,
        _disabled: { bg, borderColor, textDecoration: "none" },
      }}
      {...props}
    />
  );
};

/**
 * Inventory-styled modal with a title, body and Cancel / action footer.
 * The action button submits the form whose id is `formId`, so Enter in a field
 * submits too; pass `onAction` instead for modals without a form. Leave out
 * `actionLabel` for an information-only modal with just a Close button.
 * @param {object} props
 * @param {boolean} props.isOpen Whether the modal is shown.
 * @param {() => void} props.onClose Called on Cancel, Escape or overlay click.
 * @param {string} props.title Heading text.
 * @param {string} [props.actionLabel] Text of the action button.
 * @param {string} [props.formId] Id of the form the action button submits.
 * @param {() => void} [props.onAction] Click handler when there is no form.
 * @param {boolean} [props.isLoading] Show the action button as busy.
 * @param {boolean} [props.danger] Style the action as destructive.
 * @param {string} [props.size] Chakra modal size, "md" by default.
 */
export const InvModal = ({
  isOpen,
  onClose,
  title,
  actionLabel,
  formId,
  onAction,
  isLoading = false,
  danger = false,
  size = "md",
  children,
}) => (
  <Modal
    isOpen={isOpen}
    onClose={onClose}
    size={size}
    isCentered
    scrollBehavior="inside"
  >
    <ModalOverlay />
    <ModalContent borderRadius="12px" color={INV.ink} mx="16px">
      <ModalHeader fontSize="18px" fontWeight={600} pr="56px">
        {title}
      </ModalHeader>
      <ModalCloseButton
        w="44px"
        h="44px"
        top="8px"
        right="8px"
        _focusVisible={focusRing}
      />
      <ModalBody fontSize="14px">{children}</ModalBody>
      <ModalFooter gap="12px">
        <InvButton onClick={onClose} isDisabled={isLoading}>
          {actionLabel ? "Cancel" : "Close"}
        </InvButton>
        {actionLabel && (
          <InvButton
            primary={!danger}
            danger={danger}
            type={formId ? "submit" : "button"}
            form={formId}
            onClick={formId ? undefined : onAction}
            isLoading={isLoading}
          >
            {actionLabel}
          </InvButton>
        )}
      </ModalFooter>
    </ModalContent>
  </Modal>
);

/**
 * Page frame shared by the inventory pages: background, centred column,
 * breadcrumb, heading and header actions.
 * @param {object} props
 * @param {React.ReactNode} props.breadcrumb Small text above the heading.
 * @param {string} props.title Page heading.
 * @param {React.ReactNode} [props.actions] Buttons shown beside the heading.
 */
export const InvPage = ({ breadcrumb, title, actions, children }) => (
  <Box bg={INV.bg} color={INV.ink} minH="100vh">
    <Flex
      direction="column"
      gap="24px"
      maxW="1280px"
      mx="auto"
      px="40px"
      py="32px"
      sx={{ [MOBILE]: { px: "16px", py: "20px" } }}
    >
      <Flex
        as="header"
        justify="space-between"
        align="flex-end"
        wrap="wrap"
        gap="16px"
      >
        <Box>
          <Text m={0} mb="4px" fontSize="12px" color={INV.muted}>
            {breadcrumb}
          </Text>
          <Heading
            as="h1"
            m={0}
            fontSize="32px"
            fontWeight={600}
            lineHeight="inherit"
          >
            {title}
          </Heading>
        </Box>
        {actions && <HStack spacing="12px">{actions}</HStack>}
      </Flex>
      {children}
    </Flex>
  </Box>
);

/** Borderless accent-coloured text button, underlined on hover, focus and press. */
export const LinkButton = (props) => (
  <Button
    {...bareButton}
    px="10px"
    fontSize="13px"
    fontWeight={500}
    color={INV.accent}
    _hover={underline}
    _active={underline}
    {...props}
  />
);

/**
 * 44px square icon-only button.
 * @param {object} props
 * @param {string} props["aria-label"] Accessible name, required by IconButton.
 */
export const BareIconButton = (props) => (
  <IconButton
    {...bareButton}
    alignItems="center"
    justifyContent="center"
    w="44px"
    h="44px"
    minW="44px"
    color={INV.ink}
    {...props}
  />
);

/** White panel with a light border and rounded corners. */
export const Card = (props) => (
  <Box
    bg="white"
    border="1px solid"
    borderColor={INV.line}
    borderRadius="12px"
    {...props}
  />
);

/** Small muted secondary text, used for CAS numbers and labels. */
export const Muted = (props) => (
  <Text as="div" color={INV.muted} fontSize="12px" {...props} />
);

/**
 * Error box shown at the top of a modal when its action fails. Announced to
 * screen readers when it appears.
 */
export const ErrorMessage = (props) => (
  <Box
    role="alert"
    mb="16px"
    px="12px"
    py="10px"
    borderRadius="8px"
    bg={TONES.exp.bg}
    color={TONES.exp.color}
    fontSize="14px"
    {...props}
  />
);

/** Centred placeholder shown when a table has no rows. */
export const EmptyState = (props) => (
  <Text
    m={0}
    px="14px"
    py="32px"
    textAlign="center"
    color={INV.muted}
    {...props}
  />
);

const TONES = {
  ok: { bg: "#e3f1ea", color: "#17543d" },
  low: { bg: "#fff1d6", color: "#7a4b00" },
  exp: { bg: "#fce6e1", color: "#8a2a1b" },
  muted: { bg: "#eceeed", color: "#4a5750" },
};

/**
 * Rounded status pill.
 * @param {object} props
 * @param {"ok" | "low" | "exp" | "muted"} props.tone Colour of the pill.
 */
export const Badge = ({ tone, ...props }) => (
  <ChakraBadge
    {...TONES[tone]}
    px="9px"
    py="3px"
    borderRadius="999px"
    fontSize="12px"
    fontWeight={500}
    lineHeight="inherit"
    textTransform="none"
    whiteSpace="normal"
    {...props}
  />
);

/** Outlined chip for a single H-code. */
export const HazardTag = (props) => (
  <Box
    as="span"
    display="inline-block"
    px="6px"
    py="2px"
    mr="4px"
    border="1px solid"
    borderColor={INV.lineStrong}
    borderRadius="4px"
    fontSize="11px"
    color="#3c4c44"
    {...props}
  />
);

// Chakra's Table theme sets its own fonts, line heights and padding; these
// reset them to the mockup's values. `nested` is the denser table shown
// inside an expanded compound row.
const cellBase = {
  fontFamily: "inherit",
  lineHeight: "inherit",
  textAlign: "left",
};

/**
 * Table header cell.
 * @param {object} props
 * @param {boolean} [props.nested] Use the smaller nested-table style.
 */
export const InvTh = ({ nested = false, ...props }) => (
  <Th
    {...cellBase}
    fontWeight={600}
    fontSize={nested ? "11px" : "12px"}
    letterSpacing="0.04em"
    textTransform="uppercase"
    color={INV.muted}
    px={nested ? "12px" : "14px"}
    py={nested ? "10px" : "12px"}
    borderBottom="1px solid"
    borderColor={INV.line}
    {...props}
  />
);

/**
 * Table body cell.
 * @param {object} props
 * @param {boolean} [props.nested] Use the tighter nested-table style.
 */
export const InvTd = ({ nested = false, ...props }) => (
  <Td
    {...cellBase}
    fontSize="inherit"
    px={nested ? "12px" : "14px"}
    py={nested ? "4px" : "12px"}
    borderBottom="1px solid"
    borderColor={nested ? INV.subLine : INV.rowLine}
    verticalAlign="middle"
    {...props}
  />
);

/** Props for a Chakra Table so it matches the mockup's table base styles. */
export const tableProps = {
  variant: "unstyled",
  w: "100%",
  fontSize: "14px",
  sx: { fontVariantNumeric: "normal" },
};
