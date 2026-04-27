import { Box, Circle, Icon, Text } from "@chakra-ui/react";
import { FaBell } from "react-icons/fa";

const SIZE_STYLE_MAP = {
  xxs: { size: "16px", circlePosition: "-1.6" },
  xs: { size: "18px", circlePosition: "-2" },
  sm: { size: "20px", circlePosition: "-3" },
  md: { size: "24px", circlePosition: "-4" },
  lg: { size: "28px", circlePosition: "-6.5", textSize: "sm" },
  xl: { size: "32px", circlePosition: "-12", textSize: "lg" },
};

export const NotificationBadge = ({
  count,
  counterBg = "red.500",
  iconBg = "gray.600",
  icon = FaBell,
  size = "xs",
  ...p
}) => {
  const { size: iconSize, circlePosition } = SIZE_STYLE_MAP[size];

  return (
    <Box pos="relative" {...p}>
      <Icon as={icon} color={iconBg} boxSize={iconSize} />
      {count > 0 && (
        <Circle
          size={iconSize}
          bg={counterBg}
          color="white"
          position="absolute"
          top={circlePosition}
          right={circlePosition}
          borderColor="white"
          borderWidth={2}
        >
          <Text fontSize={SIZE_STYLE_MAP[size]?.textSize || "xs"}>{count}</Text>
        </Circle>
      )}
    </Box>
  );
};
