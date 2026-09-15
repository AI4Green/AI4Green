// src/components/Navbar.tsx
import { Flex, Text } from "@chakra-ui/react";
import { Link } from "react-router-dom";

export default function Navbar() {
  // const location = useLocation(); // to highlight the current page

  return (
    <Flex
      as="nav"
      bg="gray.500"
      color="white"
      padding={4}
      align="center"
      justify="space-between"
      boxShadow="md"
    >
      <Text fontWeight="bold" fontSize="xl">
        MyApp
      </Text>

      <Flex gap={6}>
        <nav>
          <Link to="/">Home</Link> | <Link to="/coshh">COSHH</Link> |{" "}
          {/*<Link to="/contact">Contact</Link>*/}
        </nav>
      </Flex>
    </Flex>
  );
}
