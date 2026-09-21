import {
  Button,
  HStack,
  IconButton,
  useMultiStyleConfig,
  useTab,
} from "@chakra-ui/react";
import { forwardRef } from "react";
import { FaRegTimesCircle } from "react-icons/fa";

/**
 * Custom tab component for tab list that includes a close button, which calls the onRemove function.
 * Note: This component is to be used with Chakra UI Tabs component.
 * https://v1.chakra-ui.com/docs/components/disclosure/tabs#creating-custom-tab-components
 * Props:
 * - onRemove: Function to call when the tab is removed.
 * - isRemoveDisabled: Boolean to hide the remove button.
 */
export const RemovableTab = forwardRef(function RemovableTab(
  { onRemove, isRemoveDisabled, children, ...p },
  ref,
) {
  const tabProps = useTab({ ...p, ref });
  const isSelected = !!tabProps["aria-selected"];
  const styles = useMultiStyleConfig("Tabs", tabProps);

  return (
    <HStack __css={styles.tab} {...tabProps}>
      <Button __css={styles.tab}>{children}</Button>

      {isSelected && !isRemoveDisabled && (
        <IconButton
          aria-label="Remove tab"
          isDisabled={isRemoveDisabled}
          icon={<FaRegTimesCircle />}
          variant="ghost"
          colorScheme="red"
          size="sm"
          onClick={(e) => {
            e.stopPropagation();
            onRemove();
          }}
          ml={2}
        />
      )}
    </HStack>
  );
});
