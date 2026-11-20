import { BareIconButton, INV } from "components/inventory/ui";

// smaller than inventory's 44px so step headers stay compact; still above
// WCAG 2.5.8's 24px minimum target size
const SIZE = "32px";

/**
 * Icon-only control for step cards (expand, reorder, delete, remove a
 * material): transparent until hovered, in the inventory pages' style.
 * @param {{ "aria-label": string, danger?: boolean }} props `danger` colours it
 *   for destructive actions.
 */
export const StepIconButton = ({ danger = false, ...props }) => (
  <BareIconButton
    w={SIZE}
    h={SIZE}
    minW={SIZE}
    minH={SIZE}
    borderRadius="6px"
    fontSize="13px"
    color={danger ? INV.danger : INV.ink}
    _hover={{
      bg: danger ? INV.dangerSoft : INV.accentSoft,
      _disabled: { bg: "transparent" },
    }}
    _disabled={{ opacity: 0.35, cursor: "default" }}
    {...props}
  />
);
