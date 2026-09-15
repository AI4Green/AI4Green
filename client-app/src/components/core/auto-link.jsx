import { Link } from "@chakra-ui/react";
import { Link as RouterLink } from "react-router-dom";

const isLocalUrl = (url) => url.startsWith("/");

/**
 * Figures out from the URL if the link is external or not,
 * and automatically uses a RouterLink if local
 */
export const AutoLink = ({ children, url }) => {
  const linkProps = isLocalUrl(url)
    ? { to: url, as: RouterLink }
    : { href: url, isExternal: true };
  return <Link {...linkProps}>{children}</Link>;
};
