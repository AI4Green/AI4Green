import {
  Breadcrumb,
  BreadcrumbItem,
  BreadcrumbLink,
  Text,
} from "@chakra-ui/react";
import { FaChevronRight } from "react-icons/fa";
import { Link } from "react-router-dom";

export const Breadcrumbs = ({ items = [] }) => (
  <Breadcrumb
    separator={<FaChevronRight color="gray" />}
    spacing="8px"
    padding="8px 0px"
  >
    {items.map((item, index) => (
      <BreadcrumbItem key={index} isCurrentPage={index === items.length - 1}>
        {!item.href ? (
          <Text color="gray.600">{truncateText(item?.label, 50)}</Text>
        ) : item.external ? (
          <a href={item.href} style={{ color: "#3182ce" }}>
            {truncateText(item?.label, 50)}
          </a>
        ) : (
          <BreadcrumbLink
            as={Link}
            to={item.href}
            color="blue.500"
            _hover={{ textDecoration: "underline", color: "blue.600" }}
          >
            {truncateText(item?.label, 50)}
          </BreadcrumbLink>
        )}
      </BreadcrumbItem>
    ))}
  </Breadcrumb>
);

const truncateText = (text, maxLength) => {
  if (text?.length > maxLength) {
    return text.substring(0, maxLength - 3) + "...";
  }
  return text;
};
