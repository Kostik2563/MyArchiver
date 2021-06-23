package com.shpp.p2p.cs.kscherban.assignment15;

import java.io.*;
import java.util.*;

public class Assignment15Part1 implements Constants {
    static TreeMap<Integer, String> binaryCodes;
    static File inputFile, outputFile; //input / output file
    static StringBuilder structure;

    static DataInputStream inputStream;
    static DataOutputStream outputStream;

    public static void main(String[] args) throws IOException {
        int action = determineFormatFiles(args);

        structure = new StringBuilder();
        binaryCodes = new TreeMap<>();

        inputStream = new DataInputStream(new FileInputStream(inputFile));
        outputStream = new DataOutputStream(new FileOutputStream(outputFile));
        //if archive
        if (action == 0) {
            System.out.println("Start archiving:)");
            System.out.println("Input file: " + inputFile + "\n" + "Output file: " + outputFile);
            archive();
        }
        //if unpacked
        else if (action == 1) {
            System.out.println("Start unpacking:)");
            System.out.println("Input file: " + inputFile + "\n" + "Output file: " + outputFile);
            unarchive();
        }
    }

    /**
     * determine the format of input files to know what to do archiving / unpacking
     * <p>
     * action - 0-archive / 1-unzip
     * builderInput - input file
     * builderOutput - source file
     *
     * @param args arguments
     * @return action
     */
    private static int determineFormatFiles(String[] args) {
        int action = 0;
        StringBuilder builderInput = new StringBuilder("test.txt");
        StringBuilder builderOutput = new StringBuilder(builderInput + ".par");
        if (args.length == 0)
            return action;
        String[] tmp = args[0].split("\\.");
        String last = tmp[tmp.length - 1];
        if (last.equals(PAR) || last.equals("-u")) {
            action = 1;
            if (last.equals(PAR)) {
                builderInput = new StringBuilder(args[0]);
                builderOutput = new StringBuilder(builderInput.substring(0, builderInput.length() - 4));
                if (args.length == 2)
                    builderOutput = new StringBuilder(args[1]);
            } else {
                if (args.length == 2)
                    builderInput = new StringBuilder(args[1]);
                else if (args.length == 3) {
                    builderInput = new StringBuilder(args[1]);
                    builderOutput = new StringBuilder(args[2]);
                }
            }
        } else {
            if (last.equals("-a")) {
                if (args.length == 2)
                    builderInput = new StringBuilder(args[1]);
                else if (args.length == 3) {
                    builderInput = new StringBuilder(args[1]);
                    builderOutput = new StringBuilder(args[2]);
                }
            } else {
                builderInput = new StringBuilder(args[0]);
                builderOutput = new StringBuilder(builderInput + ".par");
                if (args.length == 2) {
                    builderOutput = new StringBuilder(args[1]);
                }
            }
        }
        inputFile = new File(builderInput.toString());
        outputFile = new File(builderOutput.toString());
        return action;
    }

    /**
     * method archive the file
     *
     * @throws IOException if something went wrong :((
     * @variable time - calculation of program execution time
     * @variable lengthTree - length tree
     * @variable occurrences - symbol-frequency collection
     * @variable treeNode - collection with links to tree leaves
     * @variable tree - finally tree which is based on the Huffman algorithm
     * @variable buffer - array-buffer for byte streaming
     */
    private static void archive() throws IOException {
        long time = System.currentTimeMillis();
        int lengthTree = 0;
        HashMap<Byte, Integer> occurrences = countNumberOccurrences();
        PriorityQueue<Node> treeNode = createTree(occurrences);

        Node tree = huffman(treeNode);
        lengthTree += occurrences.size() * Byte.SIZE;

        getStructure(tree);

        lengthTree += structure.length(); //add length structure
        lengthTree += structure.length() % Byte.SIZE; //add to the length of the tree the number of "extra"
        // zeros of the last byte
        outputStream.writeShort(lengthTree);//write size tree in output file
        outputStream.write(structure.length() % Byte.SIZE); //the number of extra 0 of the last byte is written
        outputStream.writeShort(structure.length());// write length structure tree

        writeStructureToFile();

        getSymbol(tree);
        //write a collection of symbol-values
        for (Byte c : occurrences.keySet())
            binaryCodes.put(Integer.valueOf(c), tree.getBinaryCode((int) c, ""));
        inputStream = new DataInputStream(new FileInputStream(inputFile));

        StringBuilder buildByte = new StringBuilder(); //required to build bytes
        int realQuantityBytes; //the number of bytes read
        byte[] buffer = new byte[NUMBER_BYTE];

        System.out.println("******");
        for (Integer integer : binaryCodes.keySet()) {
            System.out.println(binaryCodes.get(integer));
        }
        System.out.println("******");

        while ((realQuantityBytes = inputStream.read(buffer)) != -1) {
            ArrayList<Byte> bytes = new ArrayList<>(); //buffer-bytes
            int index = 0; // index "buffer" array
            for (int i = 0; i < realQuantityBytes; i++) {
                buildByte.append(binaryCodes.get(Integer.parseInt(String.valueOf(buffer[i]))));

                while (buildByte.length() >= Byte.SIZE) {
                    bytes.add((byte) Integer.parseInt(buildByte.substring(0, Byte.SIZE), 2));
                    index++;
                    buildByte.delete(0, Byte.SIZE);
                }
            }

            //we translate a collection of bytes into an array of bytes
            byte[] bufferByteInArray = new byte[bytes.size()];
            for (int i = 0; i < bytes.size(); i++) {
                bufferByteInArray[i] = bytes.get(i);
            }

            outputStream.write(bufferByteInArray, 0, index);
        }
        int extraBitsInLastByte = 0; //extra bits (0) in last byte
        //if there are no recorded bits left
        if (buildByte.length() > 0) {
            extraBitsInLastByte = Byte.SIZE - buildByte.length();
            outputStream.write(Integer.parseInt(buildByte.toString(), 2));
        }
        outputStream.write(extraBitsInLastByte);
        System.out.println("Input file size: " + inputFile.length() + " bytes");
        System.out.println("Source file size: " + outputFile.length() + " bytes");
        System.out.println("Efficiency of archiving " + (ONE_HUNDRED - ((float) outputFile.length() /
                ((float) inputFile.length() / ONE_HUNDRED))) + "%");
        System.out.println("Execution time: " + (System.currentTimeMillis() - time) + " milliseconds");
        System.out.println("Good day!)");
    }

    /**
     * method unpacks the file
     *
     * @throws IOException if something went wrong :((
     * @variable time - calculation of program execution time
     * @variable sizeTree - the first two bytes, the size of the tree in bits
     * @variable numberBitsInLastByte - the number of bits in the last byte of the tree structure;
     * needed to restore the structure without extra zeros
     * @variable sizeStructure - the size of the tree structure without encrypted data
     * @variable keysValue - collection with a symbol and its code in a binary tree
     * @variable structure - tree structure in the form of a boolean collection
     * @variable counter - the number of bits for reading tree leaves
     * @variable symbols - tree leaves
     * @variable symbolsCopy - copy "symbols"; needed so that we can walk through the symbols and assign
     * to each of them unique binary code from the tree
     * @variable treeNode - finally tree
     * @variable extraBitsLastByte - extra 0-s in the last byte
     */
    private static void unarchive() throws IOException {
        long time = System.currentTimeMillis();
        short sizeTree = inputStream.readShort();
        int numberBitsInLastByte = inputStream.read();
        int sizeStructure = inputStream.readShort();
        HashMap<String, Byte> keysValue = new HashMap<>();

        LinkedList<Boolean> structured = getStructureTree(numberBitsInLastByte, sizeStructure);
        int counter = sizeTree - numberBitsInLastByte - sizeStructure;
        LinkedList<Integer> symbols = getSymbolTree(counter);
        LinkedList<Integer> symbolsCopy = new LinkedList<>(symbols);

        Node treeNode = restoreTree(structured, symbols);
        //write down the "key-value" collection
        for (Integer integer : symbolsCopy)
            keysValue.put(treeNode.getBinaryCode(integer, ""), (byte) Integer.parseInt(String.valueOf(integer)));
        //if there are no characters in the collection, then the file is empty and it makes no sense to run the archiver
        if (keysValue.size() == 0)
            System.exit(0);

        int extraBitsInLastByte = 0;
        int realQuantityBytes;
        StringBuilder builderByte = new StringBuilder();
        Node node = treeNode;
        byte[] buffer = new byte[NUMBER_BYTE];
        //read an array of bytes
        while ((realQuantityBytes = inputStream.read(buffer)) != -1) {
            ArrayList<Byte> bufferOutput = new ArrayList<>();
            for (int i = 0; i < realQuantityBytes; i++) {
                StringBuilder currentByte = new StringBuilder(Integer.toBinaryString(0xFF & buffer[i]));
                //if you read the last byte of encrypted data
                if (realQuantityBytes < NUMBER_BYTE && i == realQuantityBytes - 2) {
                    extraBitsInLastByte = buffer[i + 1];
                    builderByte = new StringBuilder(currentByte);
                    break;
                }

                currentByte.reverse();
                while (currentByte.length() < Byte.SIZE)
                    currentByte.append(0);
                currentByte.reverse();
                builderByte.append(currentByte); //byte in the form of a binary representation

                node = treeNode;
                int countBitsForDelete = 0; //the counter of bits to be deleted
                int endIndexForDelete = 0; //the number of the final bit to delete
                //write characters to the byte buffer
                for (int j = 0; j < builderByte.length(); j++) {
                    countBitsForDelete++;
                    node = builderByte.charAt(j) == '0' ? node.left : node.right;
                    if (node.symbol != null) {
                        bufferOutput.add((byte) Integer.parseInt(String.valueOf(node.symbol)));
                        endIndexForDelete += countBitsForDelete;
                        countBitsForDelete = 0;
                        node = treeNode;
                    }
                }
                builderByte.delete(0, endIndexForDelete);//delete bits that have already been processed
            }
            //we translate a collection of bytes into an array of bytes
            byte[] arraysByteBuffer = new byte[bufferOutput.size()];
            for (int i = 0; i < bufferOutput.size(); i++)
                arraysByteBuffer[i] = bufferOutput.get(i);

            outputStream.write(arraysByteBuffer);
        }
        //write the value of the last read byte
        builderByte.reverse();
        while (builderByte.length() < Byte.SIZE - extraBitsInLastByte)
            builderByte.append(0);
        builderByte.reverse();
        //process the last characters of the source file
        for (int j = 0; j < builderByte.length(); j++) {
            node = builderByte.charAt(j) == '0' ? node.left : node.right;
            if (node.symbol != null) {
                outputStream.write((byte) Integer.parseInt(String.valueOf(node.symbol)));
                node = treeNode;
            }
        }
        System.out.println("Execution time: " + (System.currentTimeMillis() - time) + " milliseconds");
        System.out.println("Good day!)");
    }

    /**
     * method recursively restores tree using encrypted structure, and elements of the encrypted tree
     * if true - Node is created with the children of this Node
     * if false - Node with the symbol is created
     *
     * @param structured structure in the form true/false
     * @param symbols    tree leaves
     * @return restored tree
     */
    private static Node restoreTree(LinkedList<Boolean> structured, LinkedList<Integer> symbols) {
        if (Boolean.TRUE.equals(structured.pollFirst()))
            return new Node(restoreTree(structured, symbols), restoreTree(structured, symbols));
        else
            return new Node(symbols.pollFirst());
    }

    /**
     * we read the symbol and write it in the collection
     *
     * @param counter the number of bits that are symbols; this parameter is needed to separate
     *                the encoded elements from the data itself
     * @return collection with symbols that will later be written to the tree
     * @throws IOException if something went wrong :((
     * @variable symbols - collection with elements
     */
    private static LinkedList<Integer> getSymbolTree(int counter) throws IOException {
        LinkedList<Integer> symbols = new LinkedList<>();
        while (counter > 0) {
            symbols.add(inputStream.read());
            counter -= Byte.SIZE;
        }
        return symbols;
    }

    /**
     * read the symbol from the file (! except the last) and translate it into binary form;
     * if the length is less than 8, then add the leading 0;
     * read the last byte, if the length of its binary form is less than the value of the variable passed
     * in the first parameter, then add leading zeros until the length is equal to that variable
     *
     * @param numberBitsInLastByte the number of bits in the last byte
     * @param sizeStructure        the size of the tree structure
     * @return collection with a tree structure in the form true/false
     * @throws IOException if something went wrong :((
     * @variable stringStructure - boolean collection
     * @variable structure - here we write the structure in the presented strings
     * @variable currentByteWithStructure - the current readable character that is responsible for the structure
     * @variable lastByte - the last byte read
     */
    private static LinkedList<Boolean> getStructureTree(int numberBitsInLastByte, int sizeStructure) throws IOException {
        LinkedList<Boolean> stringsStructure;
        StringBuilder structure = new StringBuilder();

        for (int i = 0; i < Math.ceil((double) sizeStructure / Byte.SIZE) - 1; i++) {
            StringBuilder currentByteWithStructure = new StringBuilder(Integer.toBinaryString(inputStream.read()));
            currentByteWithStructure.reverse();
            while (currentByteWithStructure.length() < Byte.SIZE)
                currentByteWithStructure.append(0);
            currentByteWithStructure.reverse();
            structure.append(currentByteWithStructure);
        }

        StringBuilder lastByte = new StringBuilder(Integer.toBinaryString(inputStream.read()));
        if (lastByte.length() < numberBitsInLastByte) {
            lastByte.reverse();
            while (lastByte.length() != numberBitsInLastByte)
                lastByte.append(0);
            lastByte.reverse();
        }
        structure.append(lastByte);

        stringsStructure = new LinkedList<>(translateToBoolean(structure));

        return stringsStructure;
    }

    /**
     * symbolically pass through the structure and check: if 1 is true, 0 is false
     *
     * @param structure structure tree
     * @return the structure of the tree in the form of a Boolean collection
     * @variable stringsStructure - boolean collection
     * @variable chars - structure tree in view character array
     */
    private static LinkedList<Boolean> translateToBoolean(StringBuilder structure) {
        LinkedList<Boolean> stringsStructure = new LinkedList<>();
        char[] chars = structure.toString().toCharArray();
        for (char aChar : chars) {
            if (aChar == '1')
                stringsStructure.add(true);
            else
                stringsStructure.add(false);
        }
        return stringsStructure;
    }

    /**
     * the method writes the tree structure to a file; take the first 8 characters of variable "structure",
     * form a byte and write to a file; if the characters are less than 8, then add 0 to the leading
     *
     * @throws IOException if something went wrong :((
     * @variable builderByte - stores 8 characters of variable "structure"
     */
    private static void writeStructureToFile() throws IOException {
        StringBuilder builderByte;
        while (structure.length() > Byte.SIZE) {
            builderByte = new StringBuilder(structure.substring(0, Byte.SIZE));
            outputStream.write(Integer.parseInt(builderByte.toString(), 2));
            structure.delete(0, Byte.SIZE);
        }
        if (structure.length() > 0) {
            structure.reverse();
            while (structure.length() != Byte.SIZE)
                structure.append(0);
            structure.reverse();
            outputStream.write(Integer.parseInt(structure.toString(), 2));
        }
    }

    /**
     * tree leaves are created
     *
     * @param occurrences collection with symbols and their frequency
     * @return treeNode
     * @variable treeNode - collection with tree leaves
     */
    private static PriorityQueue<Node> createTree(HashMap<Byte, Integer> occurrences) {
        PriorityQueue<Node> treeNode = new PriorityQueue<>();
        //if there are no characters, it means that the input file is empty, then an empty node is created
        if (occurrences.size() == 0)
            treeNode.add(new Node());
        for (Byte i : occurrences.keySet()) {
            treeNode.add(new Node((int) i, occurrences.get(i)));
        }
        System.out.println(treeNode.size());
        return treeNode;
    }

    /**
     * a class that creates tree leaves / knots
     *
     * @variable symbol - a character represented as a number because not all characters are printed
     * @variable frequency - frequency in the text of this symbol
     * @variable left - link to the left side of the tree
     * @variable right - link to the fight side of the tree
     */
    static class Node implements Comparable<Node> {
        Integer symbol;
        int frequency;
        Node left;
        Node right;

        public Node(Integer symbol, int frequency, Node left, Node right) {
            this.symbol = symbol;
            this.frequency = frequency;
            this.left = left;
            this.right = right;
        }

        public Node(Integer symbol, int frequency) {
            this.symbol = symbol;
            this.frequency = frequency;
        }

        public Node(Integer symbol) {
            this.symbol = symbol;
        }

        public Node(Node left, Node right) {
            this.left = left;
            this.right = right;
        }

        public Node() {
        }

        @Override
        public int compareTo(Node o) {
            return frequency - o.frequency;
        }

        // recursively we obtain a binary form (path to the symbol) of the symbol
        public String getBinaryCode(Integer character, String result) {
            if (symbol == character) {
                return result;
            } else {
                if (left != null) {
                    String path = left.getBinaryCode(character, result + 0);
                    if (path != null)
                        return path;
                }
                if (right != null)
                    return right.getBinaryCode(character, result + 1);
            }
            return null;
        }
    }

    /**
     * using the Huffman algorithm we build a tree; we take the last two elements of the sorted collection with
     * leaves and combine them with an intermediate node; repeat these steps until the size of the tree is larger than 1
     *
     * @param nodes collection with tree leaves
     * @return final tree
     */
    @SuppressWarnings("ConstantConditions")
    private static Node huffman(PriorityQueue<Node> nodes) {
        // if the leaves of the tree are less than 2, then create one empty intermediate
        // node and combine them with another intermediate node
        if (nodes.size() < 2) {
            Node first = nodes.poll();
            Node intermediate = new Node(null, first.frequency, first, new Node());
            nodes.add(intermediate);
        }
        while (nodes.size() > 1) {
            //Collections.sort(nodes);
            Node first = nodes.poll();
            Node second = nodes.poll();
            Node intermediate = new Node(null, first.frequency + second.frequency, first, second);
            nodes.add(intermediate);
        }
        return nodes.poll();
    }

    /**
     * we recursively traverse the tree by a direct algorithm and write structure tree in variable
     *
     * @param node tree
     */
    public static void getStructure(Node node) {
        if (node == null)
            return;
        /* first print data of node */
        if (node.symbol == null)
            structure.append(1);
        else
            structure.append(0);
        /* then recur on left subtree */
        getStructure(node.left);
        /* now recur on right subtree */
        getStructure(node.right);
    }

    /**
     * recursively traverse the tree in order and write to the file the leaves of the tree
     *
     * @param node tree
     * @throws IOException if something went wrong :((
     */
    private static void getSymbol(Node node) throws IOException {
        if (node == null)
            return;
        /* first recur on left child */
        getSymbol(node.left);
        /* then print the data of node */
        if (node.symbol != null) {
            outputStream.write(node.symbol);
        }
        /* now recur on right child */
        getSymbol(node.right);
    }

    /**
     * read part of the file in an array of bytes; by means of a cycle we pass on this array and we check:
     * if the collection already has such symbol we implement its value, otherwise, we add this collection
     * to a collection and we assign to it value 1
     *
     * @return method returns collection <symbol, frequency>
     * @throws IOException if something went wrong :((
     * @variable symbolOccurrences - byte collection and its frequency in the file
     * @variable realQuantityBytes - the number of bytes read
     * @variable buffer - an array-buffer with read bytes
     */
    private static HashMap<Byte, Integer> countNumberOccurrences() throws IOException {
        HashMap<Byte, Integer> symbolOccurrences = new HashMap<>();
        int realQuantityBytes;
        byte[] buffer = new byte[NUMBER_BYTE];
        while ((realQuantityBytes = inputStream.read(buffer)) != -1) {
            for (int i = 0; i < realQuantityBytes; i++) {
                if (symbolOccurrences.get(buffer[i]) != null)
                    symbolOccurrences.put(buffer[i], symbolOccurrences.get(buffer[i]) + 1);
                else
                    symbolOccurrences.put(buffer[i], 1);
            }
        }
        return symbolOccurrences;
    }
}
