import {
  Button,
  Divider,
  HStack,
  Icon,
  Text,
  useToast,
  IconButton,
} from "@chakra-ui/react";
import { FaPlus } from "react-icons/fa";
import { useSectionsListByProjectType } from "api";
import { Badge } from "components/core/Badge";
import { useBackendApi } from "contexts";
import { SECTION_TYPES, TITLE_ICON_COMPONENTS } from "constants";
import { useNavigate, useParams } from "react-router-dom";

export const BASE_PATH = "/project-type-management";

export const Area = () => {
  const navigate = useNavigate();
  const { projectTypeId, sectionId } = useParams();
  const toast = useToast();

  const { data: sections, mutate } =
    useSectionsListByProjectType(projectTypeId);
  const { sections: api } = useBackendApi();

  console.log("SSSSS", sections);

  const handleAddSection = async () => {
    try {
      const newSection = {
        projectTypeId: Number(projectTypeId),
        name: "New Section",
        sortOrder: sections ? sections.length + 1 : 1,
      };

      const response = await api.save(newSection);
      await mutate();

      if (response?.id) {
        navigate(
          `${BASE_PATH}/${projectTypeId}/sections/${response.id}?action=edit`,
        );
      }
    } catch (error) {
      toast({
        title: "Error creating section",
        status: "error",
        duration: 3000,
        isClosable: true,
      });
    }
  };

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
      <Badge label="Sections" colorScheme="blue" />
      <Divider orientation="vertical" height="20px" />

      <HStack spacing={2} align="center">
        {sections.map((section) => (
          <Button
            borderRadius="xl"
            key={section.id}
            justifyContent="flex-start"
            leftIcon={<Icon as={TITLE_ICON_COMPONENTS[section.name]} />}
            variant={Number(sectionId) === section.id ? "solid" : "ghost"}
            size="xs"
            onClick={() => {
              navigate(`${BASE_PATH}/${projectTypeId}/sections/${section.id}`, {
                replace: true,
              });
            }}
            _hover={{
              bg: "blue.50",
            }}
          >
            {section.name}
          </Button>
        ))}
        {/* Empty state helper text */}
        {sections?.length === 0 && (
          <Text fontSize="xs" color="gray.400" fontStyle="italic">
            No sections yet.
          </Text>
        )}

        <Button
          leftIcon={<FaPlus />}
          aria-label="Add Section"
          size="xs"
          colorScheme="blue"
          variant="ghost"
          borderRadius="full"
          onClick={handleAddSection}
          _hover={{ bg: "blue.50", transform: "scale(1.1)" }}
        >
          {" "}
          Add Section{" "}
        </Button>
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
