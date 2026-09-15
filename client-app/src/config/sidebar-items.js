import {
  PROJECT_TYPE_MANAGEMENT_PERMISSIONS,
  REGISTRATION_RULES_PERMISSIONS,
  USERMANAGEMENT_PERMISSIONS,
} from "constants";
import { TITLE_ICON_COMPONENTS } from "constants/experiment-ui";
import { FaUserCog } from "react-icons/fa";

export const sidebarItems = (t) => [
  {
    label: t("adminMenu.menuList.userManagement"),
    path: "/admin/user-management",
    icon: FaUserCog,
    permissions: Object.values(USERMANAGEMENT_PERMISSIONS),
  },
  {
    label: t("adminMenu.menuList.registrationRule"),
    path: "/registration-rule",
    icon: TITLE_ICON_COMPONENTS.RegistrationRule,
    permissions: [REGISTRATION_RULES_PERMISSIONS.ViewRegistrationRules],
  },
  {
    label: t("adminMenu.menuList.projectTypeManagement"),
    path: "project-type-management",
    icon: TITLE_ICON_COMPONENTS.ProjectType,
    permissions: [PROJECT_TYPE_MANAGEMENT_PERMISSIONS.ViewProjectTypes],
  },
];
