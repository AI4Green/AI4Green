import {
  Box,
  Center,
  Divider,
  Flex,
  HStack,
  Icon,
  Image,
  Link,
  Text,
  VStack,
} from "@chakra-ui/react";
import { Fragment } from "react";
import { FaReceipt } from "react-icons/fa";
import { LuShieldQuestion } from "react-icons/lu";
import { Link as RouterLink } from "react-router-dom";

const isLocalUrl = (url) => url.startsWith("/");

const FooterLink = ({ text, url, icon }) => {
  const linkProps = isLocalUrl(url)
    ? { to: url, as: RouterLink }
    : { href: url, isExternal: true };
  return (
    <Link {...linkProps}>
      <HStack spacing={2}>
        {icon && <Icon as={icon} w={4} h={4} color="gray.600" />}
        <Text fontSize={{ base: "xs", md: "sm" }} color="gray.600">
          {text}
        </Text>
      </HStack>
    </Link>
  );
};

const SmallFooter = ({ copyrightText, links }) => (
  <Flex justify="space-between" px={2}>
    {copyrightText && <Text>{copyrightText}</Text>}

    <HStack>
      {links.map((link, i) => (
        <Fragment key={i}>
          {i !== 0 && (
            <Center height="100%">
              <Divider orientation="vertical" />
            </Center>
          )}
          <FooterLink text={link.text} url={link.url} icon={link.icon} />
        </Fragment>
      ))}
    </HStack>
  </Flex>
);

const BigFooter = ({ copyrightText, links }) => (
  <VStack bg="gray.50" align="stretch" spacing={2}>
    <HStack
      align="center"
      p={4}
      w="100%"
      justify="space-between"
      direction={{ base: "column", md: "row" }}
      px={12}
    >
      <Image
        h="75px"
        src="/assets/UoN-light.png"
        alt="University of Nottingham logo"
      />

      <HStack justify="center" spacing={6}>
        {links.map((link, i) => (
          <Fragment key={i}>
            {i !== 0 && (
              <Box height={4} display={{ base: "none", md: "block" }}>
                <Divider
                  orientation="vertical"
                  borderColor="gray.400"
                  height="100%"
                />
              </Box>
            )}
            <FooterLink text={link.text} url={link.url} icon={link.icon} />
          </Fragment>
        ))}
      </HStack>
    </HStack>

    <Text fontSize="sm" color="gray.600" textAlign="center" py={6}>
      {copyrightText}
    </Text>
  </VStack>
);

export const Footer = ({ isSmall }) => {
  const copyrightText = `© ${new Date().getFullYear()} University of Nottingham. All rights reserved.`;

  const footerLinks = [
    {
      text: "Privacy Policy",
      url: "https://www.nottingham.ac.uk/utilities/privacy/privacy.aspx",
      icon: LuShieldQuestion,
    },
    {
      text: "Terms of Service",
      url: "https://www.nottingham.ac.uk/utilities/terms.aspx",
      icon: FaReceipt,
    },
  ];

  return isSmall ? (
    <SmallFooter links={footerLinks} copyrightText={null} />
  ) : (
    <BigFooter links={footerLinks} copyrightText={copyrightText} />
  );
};
