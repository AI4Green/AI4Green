import {
  Box,
  Divider,
  Drawer,
  DrawerBody,
  DrawerCloseButton,
  DrawerContent,
  DrawerOverlay,
  Grid,
  Heading,
  HStack,
  IconButton,
  useBreakpointValue,
  useDisclosure,
  VStack,
} from "@chakra-ui/react";
import { NavBar, SidebarButton } from "components/core/nav/index.js";
import { navbarItems } from "config/navbar-items";
import { sidebarItems } from "config/sidebar-items";
import { useUser } from "contexts";
import { useRef } from "react";
import { useTranslation } from "react-i18next";
import { FaBars } from "react-icons/fa";
import { Link } from "react-router-dom";

export const Sidebar = ({ children }) => {
  const { t } = useTranslation();
  const { user, isLoading } = useUser();
  const isFullMenu = useBreakpointValue({ base: false, xl: true });

  if (isLoading) return null;

  // filter items that the user has permission to access or where no permission is required
  const validItems = sidebarItems(t).filter(
    (item) =>
      !item.permissions ||
      item.permissions.every((x) => user?.permissions?.includes(x)),
  );

  // If there are no valid items, don't render the sidebar
  return validItems.length || (navbarItems.length > 0 && !isFullMenu) ? (
    <Grid templateRows="auto 1fr" minW="100%">
      <NavBar
        brand={
          <HStack spacing={2}>
            <SideMenuDrawer
              items={validItems}
              brand={
                <BrandLink
                  fontSize={{ base: "2xl", sm: "3xl" }}
                  fontWeight="hairline"
                />
              }
              isFullMenu={isFullMenu}
            />
            <BrandLink />
          </HStack>
        }
        workgroups={user.workgroups}
        user={user}
        isAuthenticated={!!user}
      />
      {children}
    </Grid>
  ) : (
    <Grid templateRows="auto 1fr" minW="100%">
      <NavBar brand={<BrandLink />} />
      {children}
    </Grid>
  );
};

const SideMenuDrawer = ({ items, brand, isFullMenu }) => {
  const DrawerState = useDisclosure();
  const btnRef = useRef();
  return (
    <>
      <IconButton
        icon={<FaBars size={20} />}
        size="xs"
        onClick={DrawerState.onOpen}
        variant="ghost"
        aria-label="Open sidebar menu"
        ref={btnRef}
      />
      {DrawerState.isOpen && (
        <Drawer
          size="xs"
          isOpen={DrawerState.isOpen}
          onClose={DrawerState.onClose}
          placement="left"
          initialFocusRef={btnRef}
        >
          <DrawerOverlay />
          <DrawerContent>
            <DrawerCloseButton aria-label="Close sidebar menu" />
            <DrawerBody>
              <VStack spacing={8} align="stretch" as="aside">
                <VStack spacing={4} mt={4} align="stretch">
                  {brand}
                  <Divider borderColor="brand.200" />
                </VStack>

                {items.length > 0 && (
                  <VStack spacing={4}>
                    {items.map((item, i) => (
                      <SidebarButton
                        item={item}
                        key={i}
                        onClose={DrawerState.onClose}
                      />
                    ))}
                    <Divider borderColor="brand.100" />
                  </VStack>
                )}

                {!isFullMenu && navbarItems.length > 0 && (
                  <VStack spacing={4}>
                    {navbarItems.map((item, i) => (
                      <SidebarButton
                        item={item}
                        key={i}
                        onClose={DrawerState.onClose}
                      />
                    ))}
                    <Divider borderColor="brand.100" />
                  </VStack>
                )}
              </VStack>
            </DrawerBody>
          </DrawerContent>
        </Drawer>
      )}
    </>
  );
};

const BrandLink = ({ ...p }) => {
  const { t } = useTranslation();
  return (
    <Link to="/">
      <Heading
        as="h1"
        fontSize={{ base: "xl", sm: "2xl" }}
        fontWeight="hairline"
        letterSpacing="tighter"
        color="brand.500"
        {...p}
      >
        {t("buttons.brand")}
      </Heading>
    </Link>
  );
};
