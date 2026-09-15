import {
  FaBook,
  FaCalculator,
  FaGithub,
  FaInfoCircle,
  FaLeaf,
} from "react-icons/fa";
import { LuLoaderPinwheel } from "react-icons/lu";

export const navbarItems = [
  {
    label: "About",
    path: "/about",
    icon: FaInfoCircle,
  },
  {
    label: "Green Chemistry",
    path: "/greenchemistry",
    icon: FaLeaf,
  },
  {
    label: "Sustainability Metrics",
    path: "/metrics",
    icon: FaCalculator,
  },
  {
    label: "Reaction Predictions",
    path: "/reaction-predictions",
    icon: LuLoaderPinwheel,
  },
  {
    label: "Documentation",
    path: "/documentation",
    icon: FaBook,
  },
  {
    label: "GitHub",
    path: "https://github.com/AI4Green/AI4Green4Students",
    icon: FaGithub,
  },
];
