import {
  Button,
  Divider,
  HStack,
  Icon,
  IconButton,
  Text,
  useToast,
  VStack,
} from "@chakra-ui/react";
import { useSectionsListByProjectType } from "api/section";
import { Badge } from "components/core/Badge";
import { InlineDraggableListField } from "components/core/forms";
import { BASE_PATH } from "components/project-type/canvas/area";
import { GLOBAL_PARAMETERS, STAGES, TITLE_ICON_COMPONENTS } from "constants";
import { useBackendApi } from "contexts";
import { Form, Formik } from "formik";
import { useRef, useState } from "react";
import { useTranslation } from "react-i18next";
import { FaPencilAlt, FaSave } from "react-icons/fa";
import { FiChevronsLeft, FiChevronsRight } from "react-icons/fi";
import { TbCancel } from "react-icons/tb";
import {
  useLocation,
  useNavigate,
  useParams,
  useSearchParams,
} from "react-router-dom";
import { array, object, string } from "yup";

import { Field } from "./field";

export const Section = ({ isCollapsed = false, projectType }) => {
  const [searchParams] = useSearchParams();
  const { projectTypeId, sectionTypeId, sectionId } = useParams();

  const navigate = useNavigate();
  const location = useLocation();

  const [feedback, setFeedback] = useState();
  const [isLoading, setIsLoading] = useState(false);
  const [isExpanded, setIsExpanded] = useState(!isCollapsed);

  const { data: sections, mutate } =
    useSectionsListByProjectType(projectTypeId);

  const { sections: api } = useBackendApi();

  const { t } = useTranslation();
  const toast = useToast();

  const formRef = useRef();

  const canEdit = projectType.stage === STAGES.Draft;
  const isEditing =
    canEdit &&
    searchParams.get("action") === "edit" &&
    searchParams.get("type") === "area-sections";

  const isEditingSectionFields =
    canEdit &&
    searchParams.get("action") === "edit" &&
    searchParams.get("type") === "section-fields";

  // const filteredSections = sections?.filter(
  //   (section) => section.sectionType.id === Number(sectionTypeId),
  // );
  const filteredSections = sections || [];

  const section = filteredSections?.find(
    (section) => section.id === Number(sectionId),
  );

  const handleSectionsSubmit = async ({ sections }) => {
    setIsLoading(true);
    const model = {
      projectTypeId: Number(projectTypeId),
      // sectionTypeId: Number(sectionTypeId),
      sections: sections.map((section) => ({
        id: section.id.startsWith("temp") ? null : Number(section.id),
        name: section.content,
        sortOrder: section.order,
      })),
    };

    try {
      await api.save(model);
      setIsLoading(false);

      toast({
        title: "Sections saved",
        status: "success",
        duration: GLOBAL_PARAMETERS.ToastDuration,
        isClosable: true,
        position: "top",
      });

      navigate(location.pathname, { replace: true });
      await mutate();
    } catch (e) {
      switch (e?.response?.status) {
        case 400: {
          setFeedback({
            status: "error",
            message: t("feedback.error_400"),
          });
          break;
        }
        default: {
          setFeedback({
            status: "error",
            message: t("feedback.error"),
          });
          break;
        }
      }
      toast({
        title: feedback.message,
        status: "error",
        duration: GLOBAL_PARAMETERS.ToastDuration,
        isClosable: true,
        position: "top",
      });
    } finally {
      setIsLoading(false);
    }
  };

  if (!sectionTypeId) return null;

  return (
    <>
      <VStack
        px={2}
        py={4}
        spacing={4}
        align="start"
        borderWidth={isExpanded ? 1 : 0}
        borderRadius={8}
        borderColor="purple.100"
      >
        <HStack justify="space-between" w="full">
          <Badge label="Sections" colorScheme="purple" />
          <HStack>
            {!isEditing && (
              <IconButton
                size="xs"
                variant="ghost"
                colorScheme="purple"
                onClick={() => setIsExpanded(!isExpanded)}
                icon={
                  <Icon
                    as={isExpanded ? FiChevronsLeft : FiChevronsRight}
                    fontSize="md"
                  />
                }
              />
            )}

            {isExpanded && canEdit && !isEditingSectionFields && (
              <Actions
                isLoading={isLoading}
                formRef={formRef}
                isEditing={isEditing}
              />
            )}
          </HStack>
        </HStack>
        {isExpanded && <Divider />}
        {!isEditing && isExpanded && <List sections={filteredSections} />}
        {isEditing && (
          <Formik
            enableReinitialize
            innerRef={formRef}
            initialValues={{
              sections:
                filteredSections?.map((section) => ({
                  id: section.id.toString(),
                  content: section.name,
                  order: section.sortOrder,
                })) || [],
            }}
            onSubmit={handleSectionsSubmit}
            validationSchema={validationSchema}
          >
            <Form>
              <InlineDraggableListField
                name="sections"
                addLabel="Add new section"
              />
            </Form>
          </Formik>
        )}
      </VStack>

      {section && <Field section={section} projectType={projectType} />}
    </>
  );
};

const List = ({ sections }) => {
  const { projectTypeId, sectionTypeId, sectionId } = useParams();
  const navigate = useNavigate();

  return (
    <VStack spacing={4} align="stretch">
      {sections
        .sort((a, b) => a.sortOrder - b.sortOrder)
        .map((section) => (
          <HStack key={section.id} spacing={2}>
            <Button
              borderRadius="xl"
              leftIcon={
                <HStack spacing={4}>
                  <Text fontSize="xs" fontWeight="light">
                    {section.sortOrder}.
                  </Text>
                  <Icon as={TITLE_ICON_COMPONENTS[section.sectionType.name]} />
                </HStack>
              }
              justifyContent="flex-start"
              w="full"
              variant={Number(sectionId) === section.id ? "solid" : "ghost"}
              size="xs"
              onClick={() => {
                navigate(
                  `${BASE_PATH}/${projectTypeId}/section-types/${sectionTypeId}/sections/${section.id}`,
                  {
                    replace: true,
                  },
                );
              }}
              _hover={{
                bg: "purple.50",
              }}
            >
              {section.name}
            </Button>
          </HStack>
        ))}
    </VStack>
  );
};

const Actions = ({ isLoading, formRef, isEditing }) => {
  const { projectTypeId, sectionTypeId } = useParams();
  const navigate = useNavigate();

  return (
    <HStack>
      {!isEditing && (
        <IconButton
          size="xs"
          icon={<FaPencilAlt />}
          aria-label="Edit"
          variant="ghost"
          colorScheme="blue"
          onClick={() => {
            navigate(
              // `${BASE_PATH}/${projectTypeId}/section-types/${sectionTypeId}/sections?action=edit&type=area-sections`,
              `${BASE_PATH}/${projectTypeId}/sections?action=edit&type=area-sections`,
              {
                replace: true,
              },
            );
          }}
        />
      )}
      {isEditing && (
        <HStack spacing={4}>
          <IconButton
            size="sm"
            icon={<FaSave />}
            aria-label="Save"
            variant="ghost"
            colorScheme="green"
            onClick={() => formRef.current.handleSubmit()}
            fontSize="lg"
            isLoading={isLoading}
          />
          <IconButton
            size="sm"
            fontSize="lg"
            icon={<TbCancel />}
            aria-label="Cancel"
            variant="ghost"
            colorScheme="yellow"
            onClick={() => {
              navigate(
                `${BASE_PATH}/${projectTypeId}/section-types/${sectionTypeId}/sections`,
                {
                  replace: true,
                },
              );
            }}
            isLoading={isLoading}
          />
        </HStack>
      )}
    </HStack>
  );
};

const validationSchema = object().shape({
  sections: array().of(
    object().shape({
      content: string()
        .required("Section name is required")
        .test("unique", "Duplicate section name", function (value) {
          if (!value) return true;

          const allSections = this.parent?.parent;

          if (!allSections || !Array.isArray(allSections)) return true;

          const normalized = value.trim().toLowerCase();
          const duplicateCount = allSections.filter(
            (section) => section.content?.trim().toLowerCase() === normalized,
          ).length;

          return duplicateCount <= 1;
        }),
    }),
  ),
});
