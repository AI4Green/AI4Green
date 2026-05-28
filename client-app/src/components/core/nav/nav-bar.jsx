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
} from "react-icons/bs";
// import { UpdateUserEmailModal } from "components/admin/user-management/modal/update-user-email";
import { LoadingModal } from "components/core/loading-modal";
import { navbarItems } from "config/navbar-items";
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

const NavButton = ({ icon, href, children }) => (
  <Button
    as="a"
    href={href}
    variant="ghost"
    color="white"
    leftIcon={<Icon as={icon} />}
    _hover={{ bg: "gray.700" }}
    _active={{ bg: "gray.800" }}
    size="sm"
    fontWeight="normal"
  >
    {children}
  </Button>
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
      bg="#777778FF"
      color="white"
      px={4}
      py={2}
      align="center"
      position="sticky"
      top={0}
      zIndex={1000}
      borderBottom="1px solid"
      borderColor="gray.700"
      boxShadow="sm"
    >
      {/* Mobile toggle */}
      {!fullMenu && (
        <IconButton
          icon={<BsList />}
          variant="ghost"
          color="white"
          mr={2}
          aria-label="Toggle navigation"
        />
      )}

      {/* Left nav */}
      {fullMenu && (
        <HStack spacing={1}>
          <Button
            as="a"
            href="/"
            variant="ghost"
            color="white"
            leftIcon={<Icon as={BsHouseDoorFill} />}
          >
            Home
          </Button>
          <Button
            as="a"
            href="/demo"
            variant="ghost"
            color="white"
            leftIcon={<Icon as={BsPencilSquare} />}
            _hover={{ bg: "gray.700" }}
            _active={{ bg: "gray.800" }}
            size="sm"
            fontWeight="normal"
          >
            Demo
          </Button>

          {/* Solvents dropdown */}
          <Menu>
            {({ isOpen }) => (
              <>
                <MenuButton
                  as={Button}
                  variant="ghost"
                  color="white"
                  size="sm"
                  leftIcon={<BsDropletFill />}
                  rightIcon={
                    <BsChevronDown
                      style={{
                        transition: "transform 0.2s",
                        transform: isOpen ? "rotate(180deg)" : "rotate(0deg)",
                      }}
                    />
                  }
                  _hover={{ bg: "gray.700" }}
                >
                  Solvents
                </MenuButton>

                <MenuList bg="gray.900" color="white">
                  <MenuItem
                    as="a"
                    href="/solvent_guide"
                    bg="gray.900"
                    color="white"
                    _hover={{
                      bg: "gray.700",
                      textDecoration: "underline",
                    }}
                  >
                    Solvent Guide
                  </MenuItem>
                  <MenuItem
                    as="a"
                    href="/solvent_PCA"
                    bg="gray.900"
                    color="white"
                    _hover={{
                      bg: "gray.700",
                      textDecoration: "underline",
                    }}
                  >
                    Solvent Surfer
                  </MenuItem>
                </MenuList>
              </>
            )}
          </Menu>

          <NavButton as="a" icon={BsArrowLeftRight} href="/retrosynthesis/">
            Retrosynthesis
          </NavButton>

          {/* Workgroup dropdown */}
          <Menu>
            <MenuButton
              as={Button}
              color="white"
              variant="ghost"
              size="sm"
              leftIcon={<BsPeopleFill />}
            >
              Workgroup
            </MenuButton>
            <MenuList bg="gray.900" borderColor="gray.700">
              {workgroups.map((wg) => (
                <MenuItem
                  key={wg}
                  as={Link}
                  color="white"
                  to={`/workgroup/${wg}`}
                >
                  {wg}
                </MenuItem>
              ))}
            </MenuList>
          </Menu>

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
          <Menu>
            <MenuButton
              as={Button}
              variant="ghost"
              color="white"
              leftIcon={<BsPersonCircle color="white" />}
              rightIcon={
                <BsChevronDown
                  style={{
                    transition: "transform 0.2s",
                  }}
                />
              }
            >
              {user.userName}
            </MenuButton>
            <MenuList bg="gray.900" color="white">
              <MenuItem
                bg="gray.900"
                color="white"
                _hover={{
                  bg: "gray.700",
                  textDecoration: "underline",
                }}
                icon={<BsBellFill />}
              >
                Notifications
              </MenuItem>
              <MenuItem
                bg="gray.900"
                color="white"
                _hover={{
                  bg: "gray.700",
                  textDecoration: "underline",
                }}
                icon={<BsGearFill />}
              >
                Manage Account
              </MenuItem>
              <MenuItem
                bg="gray.900"
                color="white"
                _hover={{
                  bg: "gray.700",
                  textDecoration: "underline",
                }}
                icon={<BsEyeFill />}
              >
                Accessibility
              </MenuItem>
              <MenuItem
                bg="gray.900"
                color="white"
                _hover={{
                  bg: "gray.700",
                  textDecoration: "underline",
                }}
                icon={<BsPeopleFill />}
              >
                Membership Summary
              </MenuItem>
              {user.role === "Admin" && (
                <MenuItem icon={<BsSpeedometer2 />}>Admin Dashboard</MenuItem>
              )}
              <MenuItem
                bg="red"
                color="white"
                _hover={{
                  bg: "gray.700",
                  textDecoration: "underline",
                }}
                icon={<BsBoxArrowRight />}
                color="white"
              >
                Logout
              </MenuItem>
            </MenuList>
          </Menu>
        ) : (
          <NavButton icon={BsBoxArrowInRight} to="/login">
            Log in
          </NavButton>
        )}
      </HStack>
    </Flex>
  );
};
