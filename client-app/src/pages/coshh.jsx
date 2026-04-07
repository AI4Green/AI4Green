import { Flex, Text } from "@chakra-ui/react";
import { useEffect, useState } from "react";

export default function COSHH() {
  const [currentTime, setCurrentTime] = useState(0);
  useEffect(() => {
    fetch("/coshh_setup/time")
      .then((res) => res.json())
      .then((data) => {
        setCurrentTime(data.time);
      });
  }, []);

  return (
    <Flex bg="blue.500" p={4} color="white" align="center">
      <Text fontSize="xl" fontWeight="bold">
        The current time is {new Date(currentTime * 1000).toLocaleString()}
      </Text>
    </Flex>
  );
}
