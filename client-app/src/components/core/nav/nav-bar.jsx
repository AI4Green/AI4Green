import {
  Avatar,
  Box,
  Button,
  HStack,
  Menu,
  MenuButton,
  MenuItem,
  MenuList,
  Text,
  useBreakpointValue,
  useDisclosure,
  VStack,
  Flex,
  IconButton,
  Image,
  MenuDivider,
  Spacer,
  Icon,
} from "@chakra-ui/react";
// import { HamburgerIcon } from "@chakra-ui/icons";
import {
  BsHouseDoorFill,
  BsPencilSquare,
  BsDropletFill,
  BsJournalText,
  BsSearch,
  BsArrowLeftRight,
  BsPeopleFill,
  BsInfoCircleFill,
  BsQuestionDiamondFill,
  BsPersonCircle,
  BsBellFill,
  BsGearFill,
  BsEyeFill,
  BsSpeedometer2,
  BsBoxArrowRight,
  BsBoxArrowInRight,
  BsList,
  BsChevronDown,
  BsExclamationDiamond,
} from "react-icons/bs";
// import { UpdateUserEmailModal } from "components/admin/user-management/modal/update-user-email";
import { LoadingModal } from "components/core/loading-modal";
import { navbarItems } from "config/navbar-items";
import { BRAND_COLORS, PROJECT_TYPE_MANAGEMENT_PERMISSIONS } from "constants";
import { useBackendApi, useUser } from "contexts";
import { forwardRef } from "react";
import { useTranslation } from "react-i18next";
import { FaEnvelope, FaHome, FaSignInAlt, FaSignOutAlt } from "react-icons/fa";
import { Link, useNavigate } from "react-router-dom";

const NavBarButton = forwardRef(function NavBarButton({ children, ...p }, ref) {
  const size = useBreakpointValue({ base: "xs", md: "sm" });
  return (
    <Button
      leftIcon={<FaHome />}
      ref={ref}
      size={size}
      borderRadius={8}
      variant="ghost"
      {...p}
      fontWeight="light"
    >
      {children}
    </Button>
  );
});

const BusyModal = ({ isOpen, onClose, verb }) => (
  <LoadingModal isOpen={isOpen} verb={verb} onClose={onClose} />
);

const LoggedInMenu = ({ user, onLogout }) => {
  const { t } = useTranslation();

  const {
    isOpen: isOpenChangeEmail,
    onOpen: onOpenChangeEmail,
    onClose: onCloseChangeEmail,
  } = useDisclosure();

  return (
    <Menu>
      <MenuButton
        py={4}
        as={Button}
        variant="ghost"
        aria-label={`User menu for user: ${user.fullName}`}
        _focus={{ boxShadow: "outline" }}
        _hover={{ backgroundColor: "none", color: "blue.500" }}
      >
        <Avatar name={user.fullName} size="sm" />
      </MenuButton>
      <MenuList color="gray.800" fontSize={{ base: "xs", md: "sm" }}>
        <MenuItem onClick={onOpenChangeEmail} icon={<FaEnvelope />}>
          Change Email
        </MenuItem>
        {isOpenChangeEmail && (
          <UpdateUserEmailModal
            isModalOpen={isOpenChangeEmail}
            onModalClose={onCloseChangeEmail}
            user={{ id: user.userId, ...user }}
          />
        )}

        <MenuItem fontSize="sm" isDisabled>
          <VStack fontSize="xxs" align="flex-start">
            <Text>{user.fullName}</Text>
            <Text>{user.email}</Text>
          </VStack>
        </MenuItem>

        <MenuItem
          onClick={onLogout}
          icon={<FaSignOutAlt />}
          color="red.600"
          _hover={{ backgroundColor: "red.100" }}
        >
          {t("buttons.logout")}
        </MenuItem>
      </MenuList>
    </Menu>
  );
};

const LoggedOutButtons = ({ t }) => (
  <NavBarButton leftIcon={<FaSignInAlt />} as={Link} to="/account/login">
    {t("buttons.login")}
  </NavBarButton>
);

const UserMenu = () => {
  const { user, signOut } = useUser();
  const { t } = useTranslation();
  const navigate = useNavigate();
  const {
    account: { logout },
  } = useBackendApi();

  const busyModalState = useDisclosure();

  const handleLogoutClick = async () => {
    busyModalState.onOpen();

    await logout();
    signOut();
    navigate("/account/login", {
      state: {
        toast: {
          title: t("logout.feedback.success"),
          status: "success",
          duration: 2500,
          isClosable: true,
        },
      },
    });

    busyModalState.onClose();
  };

  return (
    <>
      {user ? (
        <Box>
          <LoggedInMenu user={user} onLogout={handleLogoutClick} />
        </Box>
      ) : (
        <LoggedOutButtons t={t} />
      )}
      <BusyModal
        isOpen={busyModalState.isOpen}
        onClose={busyModalState.onClose}
        verb={t("logout.feedback.busy")}
      />
    </>
  );
};

// hover, keyboard focus and open-dropdown states share the green text and underline bar
const navItemStyles = {
  variant: "ghost",
  size: "sm",
  fontWeight: "normal",
  color: BRAND_COLORS.textMuted,
  _before: {
    position: "absolute",
    left: 0,
    bottom: "-4px",
    width: "100%",
    height: "3px",
    bg: BRAND_COLORS.green,
    borderRadius: "1.5px",
  },
  _hover: {
    bg: "transparent",
    color: BRAND_COLORS.green,
    _before: { content: '""' },
  },
  _focusVisible: {
    color: BRAND_COLORS.green,
    _before: { content: '""' },
  },
  _expanded: {
    color: BRAND_COLORS.green,
    _before: { content: '""' },
  },
  _active: { bg: "transparent" },
};

const menuListStyles = {
  bg: "white",
  color: BRAND_COLORS.navy,
  borderColor: BRAND_COLORS.border,
  borderRadius: "8px",
  boxShadow: "0 8px 24px rgba(24, 36, 44, 0.08)",
  py: 2,
  minW: "250px",
};

const menuItemStyles = {
  as: "a",
  px: 6,
  py: 2,
  bg: "white",
  color: BRAND_COLORS.navy,
  _hover: {
    bg: BRAND_COLORS.greenTint,
    color: BRAND_COLORS.green,
    fontWeight: "bold",
  },
  _focus: {
    bg: BRAND_COLORS.greenTint,
    color: BRAND_COLORS.green,
    fontWeight: "bold",
  },
};

const dangerMenuItemStyles = {
  ...menuItemStyles,
  color: BRAND_COLORS.danger,
  _hover: { bg: BRAND_COLORS.dangerTint, color: BRAND_COLORS.danger },
  _focus: { bg: BRAND_COLORS.dangerTint, color: BRAND_COLORS.danger },
};

const DropdownChevron = ({ isOpen }) => (
  <BsChevronDown
    style={{
      transition: "transform 0.2s",
      transform: isOpen ? "rotate(180deg)" : "rotate(0deg)",
    }}
  />
);

const NavButton = ({ icon, href, children }) => (
  <Button as="a" href={href} leftIcon={<Icon as={icon} />} {...navItemStyles}>
    {children}
  </Button>
);

const NavDropdown = ({ icon, label, children }) => (
  <Menu>
    {({ isOpen }) => (
      <>
        <MenuButton
          as={Button}
          leftIcon={<Icon as={icon} />}
          rightIcon={<DropdownChevron isOpen={isOpen} />}
          {...navItemStyles}
        >
          {label}
        </MenuButton>
        <MenuList {...menuListStyles}>{children}</MenuList>
      </>
    )}
  </Menu>
);

export const NavBar = ({
  brand,
  user,
  workgroups,
  notifications,
  isAuthenticated,
}) => {
  const fullMenu = useBreakpointValue({ base: false, md: true });

  return (
    <Flex
      as="nav"
      bg="white"
      color={BRAND_COLORS.textMuted}
      px={4}
      py={2}
      align="center"
      position="sticky"
      top={0}
      zIndex={1000}
      borderBottom="1px solid"
      borderColor={BRAND_COLORS.border}
      boxShadow="0 2px 8px rgba(24, 36, 44, 0.04)"
    >
      {/* Mobile toggle */}
      {!fullMenu && (
        <IconButton
          icon={<BsList />}
          variant="ghost"
          color={BRAND_COLORS.navy}
          mr={2}
          aria-label="Toggle navigation"
        />
      )}

      {/* served by Flask, like the logo in the footer */}
      <Box as="a" href="/" mr={6} flexShrink={0}>
        <Image
          h="40px"
          w="auto"
          src="/static/img/AI4Green_logo_light_small.png"
          alt="AI4Green home"
        />
      </Box>

      {/* Left nav */}
      {fullMenu && (
        <HStack spacing={1}>
          <NavButton icon={BsHouseDoorFill} href="/">
            Home
          </NavButton>
          <NavButton icon={BsPencilSquare} href="/demo">
            Demo
          </NavButton>

          <NavDropdown icon={BsDropletFill} label="Solvents">
            <MenuItem
              icon={<BsJournalText />}
              href="/solvent_guide"
              {...menuItemStyles}
            >
              Solvent Guide
            </MenuItem>
            <MenuItem
              icon={<BsSearch />}
              href="/solvent_PCA"
              {...menuItemStyles}
            >
              Solvent Surfer
            </MenuItem>
          </NavDropdown>

          <NavButton icon={BsArrowLeftRight} href="/retrosynthesis/">
            Retrosynthesis
          </NavButton>

          <NavDropdown icon={BsPeopleFill} label="Workgroup">
            {workgroups.map((wg) => (
              <MenuItem
                key={wg}
                href={`/workgroup/${encodeURIComponent(wg)}`}
                {...menuItemStyles}
              >
                {wg}
              </MenuItem>
            ))}
          </NavDropdown>

          <NavButton icon={BsSearch} href="/search">
            Search
          </NavButton>
        </HStack>
      )}

      <Spacer />

      {/* Reaction save indicator placeholder */}
      <Box id="reaction-saved-indicator" mx={4} />

      {/* Right nav */}
      <HStack spacing={1}>
        <NavButton icon={BsInfoCircleFill} href="/about">
          About AI4Green
        </NavButton>

        <NavButton icon={BsQuestionDiamondFill} href="/info">
          Help
        </NavButton>

        {isAuthenticated ? (
          <NavDropdown icon={BsPersonCircle} label={user.userName}>
            <MenuItem
              icon={<BsBellFill />}
              href="/notifications"
              {...menuItemStyles}
            >
              Notifications
            </MenuItem>
            <MenuItem
              icon={<BsGearFill />}
              href="/manage_account"
              {...menuItemStyles}
            >
              Manage Account
            </MenuItem>
            <MenuItem
              icon={<BsEyeFill />}
              href="/accessibility"
              {...menuItemStyles}
            >
              Accessibility
            </MenuItem>
            <MenuItem
              icon={<BsPeopleFill />}
              href="/workgroup_membership_summary"
              {...menuItemStyles}
            >
              Membership Summary
            </MenuItem>
            {/* COSHH pages are SPA routes, so navigate client side */}
            <MenuItem
              icon={<BsExclamationDiamond />}
              {...menuItemStyles}
              as={Link}
              to="/my-coshh-forms"
            >
              My COSHH Forms
            </MenuItem>
            {/* only H&S Admins are given CreateProjectTypes by /api/users/me */}
            {user.permissions?.includes(
              PROJECT_TYPE_MANAGEMENT_PERMISSIONS.CreateProjectTypes,
            ) && (
              <MenuItem
                icon={<BsExclamationDiamond />}
                {...menuItemStyles}
                as={Link}
                to="/coshh-form-management"
              >
                Manage COSHH Templates
              </MenuItem>
            )}
            {user.role === "Admin" && (
              <MenuItem
                icon={<BsSpeedometer2 />}
                href="/admin_dashboard"
                {...menuItemStyles}
              >
                Admin Dashboard
              </MenuItem>
            )}
            <MenuDivider borderColor={BRAND_COLORS.border} />
            <MenuItem
              icon={<BsBoxArrowRight />}
              href="/auth/logout"
              {...dangerMenuItemStyles}
            >
              Logout
            </MenuItem>
          </NavDropdown>
        ) : (
          <NavButton icon={BsBoxArrowInRight} href="/auth/login">
            Log in
          </NavButton>
        )}
      </HStack>
    </Flex>
  );
};
