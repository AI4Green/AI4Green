import {
  Badge as ChakraBadge,
  Box,
  Button,
  IconButton,
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

/**
 * Outlined button, or the solid accent button with `primary`.
 * @param {object} props
 * @param {boolean} [props.primary] Use the solid accent style.
 */
export const InvButton = ({ primary = false, ...props }) => {
  const bg = primary ? INV.accent : "white";
  const borderColor = primary ? INV.accent : INV.control;
  return (
    <Button
      {...bareButton}
      px="16px"
      borderRadius="8px"
      border="1px solid"
      borderColor={borderColor}
      bg={bg}
      color={primary ? "white" : INV.ink}
      fontSize="14px"
      fontWeight={primary ? 500 : 400}
      _disabled={{ opacity: 0.5, cursor: "default" }}
      _hover={{
        bg: primary ? INV.accentHover : INV.sub,
        borderColor: primary ? INV.accentHover : INV.accent,
        _disabled: { bg, borderColor },
      }}
      {...props}
    />
  );
};

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
