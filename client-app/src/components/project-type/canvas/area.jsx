import { Button, Divider, HStack, Icon } from "@chakra-ui/react";
import { useSectionTypesList } from "api/project-type";
import { Badge } from "components/core/Badge";
import { SECTION_TYPES, TITLE_ICON_COMPONENTS } from "constants";
import { useNavigate, useParams } from "react-router-dom";

export const BASE_PATH = "/project-type-management";

export const Area = () => {
  const navigate = useNavigate();
  const { projectTypeId, sectionTypeId } = useParams();
  const { data: sectionTypes } = useSectionTypesList();
  console.log("SECTION TYPES", sectionTypes);

  return (
    <HStack
      p={2}
      spacing={4}
      align="center"
      borderWidth={1}
      borderRadius={20}
      borderColor="blue.100"
      maxW="2xl"
    >
      <Badge label="Areas" colorScheme="blue" />

      <HStack spacing={2} align="center">
        <Divider orientation="vertical" height="20px" />
        {sectionTypes.map((sectionType) => (
          <Button
            borderRadius="xl"
            key={sectionType.id}
            justifyContent="flex-start"
            leftIcon={<Icon as={TITLE_ICON_COMPONENTS[sectionType.name]} />}
            variant={
              Number(sectionTypeId) === sectionType.id ? "solid" : "ghost"
            }
            size="xs"
            onClick={() => {
              navigate(
                `${BASE_PATH}/${projectTypeId}/section-types/${sectionType.id}/sections`,
                {
                  replace: true,
                },
              );
            }}
            _hover={{
              bg: "blue.50",
            }}
          >
            {sectionType.name}
          </Button>
        ))}
      </HStack>
    </HStack>
  );
};

const sectionTypeLabels = {
  [SECTION_TYPES.LiteratureReview]: "Literature Review",
  [SECTION_TYPES.Plan]: "Plan",
  [SECTION_TYPES.Note]: "Note",
  [SECTION_TYPES.Report]: "Report",
  [SECTION_TYPES.ProjectGroup]: "Project Group Activities",
};
