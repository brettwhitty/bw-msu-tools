#!/usr/bin/perl

use strict;
use warnings;

##
## create database solcomp_mail;
## create table mail (id serial primary key, address text, hash text, verified boolean default false);
## CREATE UNIQUE INDEX idx_address ON mail(address);
## CREATE UNIQUE INDEX idx_hash ON mail(hash);
##

use DBI;
use Email::Simple;
use Email::Simple::Creator;
use Email::Send;

use Math::Random;
use Digest::MD5 qw{ md5_hex };

my $db_user = $ENV{'DB_USER'};
my $db_pass = $ENV{'DB_PASS'};

my $dbh = DBI->connect("dbi:Pg:dbname=solcomp_mail;host=pg", $db_user, $db_pass);

my $sth;

$sth = $dbh->prepare('select * from mail where hash is null');
$sth->execute();
while (my $row = $sth->fetchrow_hashref) {
    my $hash = md5_hex( reverse($row->{'address'}) . time() );
    send_email_notice($row->{'address'}, $hash);
    my $uph = $dbh->prepare('update mail set hash = ? where address = ?');
    $uph->execute($hash, $row->{'address'});
}

sub send_email_notice {
    my ($address, $hash) = @_;

    my $body = <<BODY_TEXT;
Hello,

You are receiving this email in response to a request from your email address to sign up for the Solanaceae Genomics Resource mailing list.

To activate your subscription to the mailing list, please visit the following link:
http://solanaceae.plantbiology.msu.edu/cgi-bin/contact/mailing_list.cgi?id=$hash

If you have not subscribed on our site at:
http://solanaceae.plantbiology.msu.edu
then you are receiving this email in error. If you ignore this email, you will not receive any further emails.

Regards,

The Solanaceae Genomics Resource Team
sgr\@plantbiology.msu.edu
BODY_TEXT

    my $email = Email::Simple->create(
        header  =>  [
            From    =>  'sgr@plantbiology.msu.edu',
            To      =>  $address,
            Subject =>  'Confirm Subscription to Solanaceae Genomics Resource Mailing List',
        ],
        body    => $body,
    );
   
    my $sender = new Email::Send;
    $sender->send($email);
}
