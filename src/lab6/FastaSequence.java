package lab6;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;

public class FastaSequence
{
	private String header;
	private String sequence;

	// constructor
	public FastaSequence(String header, String sequence)
	{

		this.header = header;
		this.sequence = sequence;
	}

	// static factory method
	public static List<FastaSequence> readFastaFile(String filePath) throws Exception
	{
		// generate a list to store header/sequence
		List<FastaSequence> list = new ArrayList<FastaSequence>();
		// read fasta file
		BufferedReader reader = new BufferedReader(new FileReader(new File(filePath)));
		// check if the file is a fasta file
		// read the first line
		String firstLine = reader.readLine();
		String header;
		StringBuffer sequenceb = new StringBuffer();
		if (firstLine.startsWith(">"))
		{
			header = firstLine.substring(1).trim();

		} else
		{
			reader.close();
			throw new Exception("Please make sure this is a fasta file!");

		}
		for (String nextLine = reader.readLine(); nextLine != null; nextLine = reader.readLine())
		{
			// read from the second line
			// the program will implement else statement first
			if (nextLine.startsWith(">"))
			{
				// whenever reaches a new header, save the previous header/sequence
				FastaSequence pair = new FastaSequence(header, sequenceb.toString());
				list.add(pair);
				// update the header with new header
				header = nextLine.substring(1).trim();
				// empty the sequence stringBuffer
				sequenceb.setLength(0);
			} else
			{
				// save 1-many lines of sequence into the stringBuffer
				sequenceb.append(nextLine.trim());
			}
		}
		// generate a new FastaSequence to store the last header/sequence
		list.add(new FastaSequence(header, sequenceb.toString()));
		// close the reader
		reader.close();
		// return the list
		return list;

	}

	/*
	 * add a static method to do this writes each unique sequence to the output file
	 * with the # of times each sequence was Seen in the input file as the header
	 * (sorted with the sequence seen the fewest times the first)
	 */

	public static void writeUnique(File inFile, File outFile) throws Exception
	{
		// generate a list of FastaSequence
		List<FastaSequence> fastaList = FastaSequence.readFastaFile(inFile.toString());
		// read the list and store unique sequences and corresponding counts in HashMap
		HashMap<String, Integer> seqCount = new HashMap<String, Integer>();
		for (FastaSequence fs : fastaList)
		{
			String key = fs.getSequence();
			if (!seqCount.containsKey(key))
			{
				Integer count = 1;
				seqCount.put(key, count);
			} else
			{
				Integer count = seqCount.get(key);
				count++;
				seqCount.put(key, count);
			}
		}
		// sort by the value of the hashMap
		Set<Entry<String, Integer>> mySet = seqCount.entrySet();
		// Sort method needs a list, so convert the set to a list first
		List<Entry<String, Integer>> myList = new ArrayList<Entry<String, Integer>>(mySet);
		// define comparator
		Comparator<Entry<String, Integer>> valueComparator = new Comparator<Entry<String, Integer>>()
		{
			@Override
			public int compare(Entry<String, Integer> o1, Entry<String, Integer> o2)
			{
				// in ascending order
				return o1.getValue().compareTo(o2.getValue());
			}
		};

		// Sort HashMap by value using new defined comparator
		Collections.sort(myList, valueComparator);
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		for (Entry<String, Integer> entry : myList)
		{
			writer.write(">" + entry.getValue() + "\n");
			writer.write(entry.getKey() + "\n");
		}
		// flush the buffer
		writer.flush();
		// close the writer
		writer.close();

	}

	/*
	 * Add a static method to perform this: take that Fasta file and produce a
	 * reduced view. This will be a spreadsheet counting the number of times we’ve
	 * seen each unique sequence with the DNA sequences as column headers and the
	 * samples as rows…
	 */

	public static void writeFileInReducedView(File inFile, File outFile) throws Exception
	{
		// generate a list of FastaSequence
		List<FastaSequence> fastaList = FastaSequence.readFastaFile(inFile.toString());
		// save in a LinkedHashMap to keep the order of sequence as insertion-order
		// generate the key by merging sequence and sample
		// while the value is the number of this sequence occurs in this sample
		LinkedHashMap<String, Integer> map = new LinkedHashMap<String, Integer>();
		// store seq_sam:count into the map
		for (FastaSequence fs : fastaList)
		{
			// get sequence as part of the final key
			String seqKey = fs.getSequence();
			// get the sample ID from header by applying StringTokenizer class as the second
			// part of final key
			StringTokenizer sToken = new StringTokenizer(fs.getHeader());
			sToken.nextToken();
			String samKey = sToken.nextToken();
			// generate the final key
			String key = seqKey + "~" + samKey;
			if (!map.containsKey(key))
			{
				Integer count = 1;
				map.put(key, count);
			} else
			{
				Integer count = map.get(key);
				count++;
				map.put(key, count);
			}
		}
		// store sequence into a linkedHashSet for the column name in the output file
		LinkedHashSet<String> column = new LinkedHashSet<String>();
		// store sample into a treeSet as the row names to get the natural order
		TreeSet<String> row = new TreeSet<String>();
		for (String ckey : map.keySet())
		{
			// get sequence
			column.add(ckey.substring(0, ckey.lastIndexOf("~")));
			// get sample
			row.add(ckey.substring(ckey.lastIndexOf("~") + 1));
		}
		// write into the output
		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));
		// write the first line
		String firstLine = "sample" + "\t" + String.join("\t", column) + "\n";
		writer.write(firstLine);
		// write 2nd-last line
		// first loop through rows
		Iterator<String> itr = row.iterator();
		while (itr.hasNext())
		{
			String element = itr.next();
			writer.write(element + "\t");
			// then loop through column
			Iterator<String> itc = column.iterator();
			while (itc.hasNext())
			{
				String fkey = itc.next() + "~" + element;
				if (map.containsKey(fkey))
				{
					writer.write(map.get(fkey) + "\t");
				} else
				{
					writer.write(0 + "\t");
				}
			}
			writer.write("\n");
		}
		writer.flush();
		writer.close();
	}

	// returns the header of this sequence without the “>”
	public String getHeader()
	{
		return header;
	}

	// returns the DNA sequence of this FastaSequence
	public String getSequence()
	{
		return sequence;
	}

	// returns the number of G’s and C’s divided by the length of this sequence
	public float getGCRatio()
	{
		int countGC = 0;

		String currentSequence = this.getSequence().toUpperCase();

		for (int x = 0; x < currentSequence.length(); x++)
		{
			char target = currentSequence.charAt(x);

			if (target == 'C' || target == 'G')
				countGC++;
		}

		return (float) countGC / currentSequence.length();
	}

	public static void main(String[] args) throws Exception
	{
		// ask user for the absolute path of a fasta file
		System.out.println("Please type the absolute path of your input fasta file");
		String filePath = System.console().readLine();
		File fileIn = new File(filePath);
		// ask user for the absolute path of a fasta file
		System.out.println("Please type the absolute path of your output fasta file");
		File fileOut = new File(System.console().readLine());

		// generate FastaSequence
		List<FastaSequence> fastaList = FastaSequence.readFastaFile(filePath);
		// Print out header-sequence
		for (FastaSequence fs : fastaList)
		{
			System.out.println(fs.getHeader());
			System.out.println(fs.getSequence());
			System.out.println(fs.getGCRatio());
		}
		// write sorted count-sequence in the ascending order of count into a file
		// FastaSequence.writeUnique(fileIn, fileOut);
		// generate a reduce-view of sequence counts in each sample
		FastaSequence.writeFileInReducedView(fileIn, fileOut);

	}
}