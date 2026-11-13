import {
  Box,
  Flex,
  Heading,
  HStack,
  Image,
  Link,
  Text,
} from "@chakra-ui/react";
import { SINGULAR_ROLE_LABELS } from "components/core/structure-sketcher/structure-search-labels.js";
import { Badge, Card, INV, Muted } from "components/inventory/ui.jsx";

/**
 * Builds the link to open a reaction in the (Flask) reaction constructor.
 * Every segment is encoded, so names containing "/" or "?" can't change the
 * route.
 * @param {object} reaction A search result.
 * @returns {string} The URL.
 */
const reactionUrl = ({ workgroup, workbook, reactionId }) =>
  `/sketcher/${[workgroup, workbook, reactionId]
    .map(encodeURIComponent)
    .join("/")}/no`;

/**
 * Formats an ISO date for display, e.g. "11 Nov 2026".
 * @param {string | null} iso The date.
 * @returns {string} The formatted date, or "" if there isn't one.
 */
const formatDate = (iso) =>
  iso
    ? new Date(iso).toLocaleDateString(undefined, {
        day: "numeric",
        month: "short",
        year: "numeric",
      })
    : "";

/**
 * Whether an image source is safe to show: only image data URLs are, never
 * e.g. a javascript: URL. The backend checks this too.
 * @param {string | null} src The image source.
 * @returns {boolean} True if it's an image data URL.
 */
const isImageDataUrl = (src) =>
  typeof src === "string" && src.startsWith("data:image/");

/**
 * One reaction in the search results. All text is rendered as text, never as
 * HTML.
 * @param {object} props
 * @param {object} props.reaction A result from the search API.
 */
export const ReactionResultCard = ({ reaction }) => (
  <Card as="article" p="16px">
    <Flex gap="16px" direction={{ base: "column", md: "row" }}>
      <Box
        flexShrink={0}
        w={{ base: "100%", md: "260px" }}
        h="140px"
        bg={INV.sub}
        borderRadius="8px"
        display="flex"
        alignItems="center"
        justifyContent="center"
        overflow="hidden"
      >
        {isImageDataUrl(reaction.image) ? (
          <Image
            src={reaction.image}
            alt={`Reaction scheme for ${reaction.name}`}
            maxW="100%"
            maxH="100%"
            objectFit="contain"
          />
        ) : (
          <Muted>No reaction scheme</Muted>
        )}
      </Box>

      <Box flex={1} minW={0}>
        <Flex justify="space-between" align="flex-start" gap="12px" wrap="wrap">
          <Box minW={0}>
            <Heading as="h3" m={0} fontSize="16px" fontWeight={600}>
              {reaction.name}
            </Heading>
            <Muted>
              {reaction.reactionId} · {reaction.workgroup} / {reaction.workbook}
            </Muted>
          </Box>
          <Link
            href={reactionUrl(reaction)}
            target="_blank"
            rel="noopener noreferrer"
            color={INV.accent}
            fontWeight={500}
            fontSize="14px"
            aria-label={`Open ${reaction.name} in a new tab`}
          >
            Open reaction
          </Link>
        </Flex>

        {reaction.description && (
          <Text mt="8px" mb={0} fontSize="14px" noOfLines={3}>
            {reaction.description}
          </Text>
        )}

        <HStack mt="10px" spacing="8px" wrap="wrap">
          <Badge
            tone={reaction.completionStatus === "complete" ? "ok" : "muted"}
          >
            {reaction.completionStatus === "complete"
              ? "Complete"
              : "Not complete"}
          </Badge>
          {reaction.matchedRoles?.map((role) => (
            <Badge key={role} tone="low">
              Matched in {SINGULAR_ROLE_LABELS[role] ?? role}
            </Badge>
          ))}
          <Muted as="span">
            Created {formatDate(reaction.timeOfCreation)}
            {reaction.creatorUsername && ` by ${reaction.creatorUsername}`}
          </Muted>
        </HStack>
      </Box>
    </Flex>
  </Card>
);
